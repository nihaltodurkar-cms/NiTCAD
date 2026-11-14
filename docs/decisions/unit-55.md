<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 55 — theoretical calculation validation framework: research and design proposal (proposed 2026-10-31; on `verify/theoretical-validation`)

**Owner request (2026-10-31):** "Design a comprehensive Theoretical Calculation Validation Framework for NiTCAD.
Objective: Establish an independent theoretical reference system to verify the mathematical, physical and numerical
correctness of simulation results." Ten requirements (R1–R10 below) and five deliverables. "Initially, produce only the
research and design proposal. Do not write implementation code or modify existing architecture. Prioritise theoretical
rigour, physical validity, numerical accuracy and independent verification over simply achieving passing tests."

**Owner rulings (2026-10-31):** "Approve TVQ1–TVQ6 with the following conditions: 1. Consistency checks must support
`NOT APPLICABLE` where physical assumptions do not apply. 2. Error decomposition must account for interactions between
error sources; do not assume contributions are strictly additive. 3. Every theoretical reference must have documented
provenance, independent derivation and verification. For TVQ4, prefer an independently implemented symmetric
tridiagonal eigensolver, with Eigen used as a secondary cross-check. Start with TVR-0 (existing validation audit)."
Where each condition lands: (1) §4.4 outcomes and §5.1 (TVD-16); (2) §4.1 and §4.3 (TVD-17); (3) §3.3 and §2.2
(TVD-18); TVQ4 in TVB-30 and TVD-3. The TVR-0 audit is §10 below. "Documentation and research only. No implementation
code, architecture changes or commits until review."

**Scope of this entry.** A proposal only: no code, no change to ARCHITECTURE.md, no change to any existing test, gate or
tolerance. It plans an addition (a new `verify` layer and its reference library) that consumes results through the
public API, as ARCHITECTURE 6.12.6 already foresees ("a future verification layer … consumes results rather than
modifying the solver"). Where a benchmark needs a capability NiTCAD lacks, it is listed as a dependency (TVD-15), not
designed here.

**Requirements map (owner R1–R10 → sections):** R1 governing equations → §1; R2 independent references → §2–§3;
R3 assumptions and limits → §1 and the catalogue's columns; R4 comparison → §4; R5 error metrics → §4; R6 consistency,
conservation, convergence → §5; R7 exact vs approximate vs empirical → §1.2; R8 criteria and reproducible benchmarks →
§4, §6; R9 limitations → §7; R10 extensibility → §3.4. Deliverables: framework and methodology (§1–§2), reference
architecture (§3), benchmark catalogue (§6), metrics and acceptance (§4), roadmap (§8).

#### 1. Theoretical framework and methodology

**1.1 What is being verified.** Following the standard verification-and-validation hierarchy (Roache 1998; Oberkampf
and Roy 2010; ASME V&V 20), the framework separates five questions, which must be answered in this order, because an
answer to a later one is meaningless while an earlier one is open:

| Level | Question | Evidence | NiTCAD today |
|---|---|---|---|
| L0 Reference verification | Is the reference itself computed correctly? | The reference's own convergence; two independent references agreeing; published tabulated values | Not systematic |
| L1 Model-implementation verification | Does the code evaluate the model equations it claims (each mobility, recombination, statistics, band, tunnelling formula)? | Formula and table comparison to the cited source | Per unit, partial (unit tests against hand values) |
| L2 Code verification | Does the discretized solver converge to the exact solution of the continuous model, at the expected order? | Exact solutions, manufactured solutions, discrete exactness properties | Per unit, partial (Unit 16's H1/L2 orders, Kirchhoff U4, a few others) |
| L3 Solution verification | How large is the numerical error of *this* result (mesh, time step, Newton, linear solve, rounding)? | Convergence studies, Richardson / GCI, iterative-error estimates | Per unit; Unit 32's adjoint estimate for currents; Unit 17's proposed AD11 |
| L4 Theory comparison | How far is the model's answer from established device theory in that theory's domain of validity? | Analytical and semi-analytical device theory (Shockley, Pao–Sah, …) | Scattered (legacy MOS-C landmarks, ideal-diode law) |
| L5 Validation | How far is the model from measurement? | Published data (Unit 33) | Unit 33's baseline (ARCHITECTURE 10.1) |

This framework covers L0–L4 and catalogues L5 (Unit 33) without changing it. **Rule TVM1:** a discrepancy at level k
is not attributed to a level above k (e.g. to "model error") until levels below it pass for that benchmark.

**1.2 Classes of reference (R7).** Every reference carries exactly one class, and the class fixes what a disagreement
can mean and which criterion applies (§4):

| Class | Definition | A disagreement means | Example |
|---|---|---|---|
| **E — exact** | An exact solution of the *same continuous model* NiTCAD solves (closed form, or semi-analytic: reduced to quadrature, a scalar root, or a 1D ODE solved to a certified tolerance far below the test's) | A code error, or numerical error not yet converged | Equilibrium Poisson–Boltzmann first integral; MOS-C Q_s(ψ_s); manufactured solutions |
| **D — discrete-exact** | A property the *discrete* scheme satisfies exactly, so the result must hold to rounding | A code error | Scharfetter–Gummel exactness at constant field and current; zero current at equilibrium; Gauss's law on the dual mesh |
| **A — approximate theory** | An established approximation (depletion, low injection, gradual channel, WKB, charge sheet) with a stated domain of validity | Approximation error (expected, to be bounded), or a code error if beyond the bound | Shockley diode, depletion capacitance, subthreshold swing, Fowler–Nordheim |
| **M — empirical model** | A fitted formula with published parameters | Only an implementation or transcription error (L1); never evidence about physics | Masetti, Caughey–Thomas, Canali, Lombardi, Slotboom, Kane, van Overstraeten–de Man |
| **H — higher-fidelity theory** | A more complete physical theory than the model (Schrödinger–Poisson, exact Airy-function transmission) | Model-form error of NiTCAD's model, quantified, not a code error | MLDA vs 1D Schrödinger–Poisson; WKB vs exact transmission |
| **P — published numerical** | Another code's published numerical result | A discrepancy to investigate; weak evidence on its own | Aste–Vahldieck ballistic diode (Unit 52) |
| **X — experimental** | Measured data with uncertainty (L5) | Model or parameter error | Unit 33's resistivity and n_i |

**1.3 Governing equations in scope (R1).** Each benchmark names the subset it exercises, the boundary conditions and
the assumptions (catalogue columns). The equations, as NiTCAD solves them (ARCHITECTURE 5–6):

| # | Equation / model | Form | Where in NiTCAD | Reference classes available |
|---|---|---|---|---|
| G1 | Poisson | ∇·(ε∇ψ) = −q(p − n + N_D⁺ − N_A⁻) + interface and fixed charge | assemble (box method) | E, D |
| G2 | Carrier continuity (steady, transient) | ∂n/∂t = (1/q)∇·J_n + G − R; ∂p/∂t = −(1/q)∇·J_p + G − R | assemble (SG fluxes) | E (MMS, 1D linearised), D, A |
| G3 | Drift-diffusion currents | J_n = qμ_n n E + qD_n∇n (Einstein relation; FD-generalised) | SG flux | D, E |
| G4 | Statistics | Boltzmann; Fermi–Dirac (F_{1/2}), incomplete ionization | physics | E (integrals), M |
| G5 | Band structure | E_g(T), χ, N_c, N_v (Couderc, legacy), BGN (Slotboom) | physics | M |
| G6 | Mobility | Caughey–Thomas, Masetti, Canali v(E), Lombardi | physics | M (+ E for uniform-field drift) |
| G7 | Recombination | SRH, Auger, radiative, surface, interface traps (steady and dynamic) | physics, assemble | M, E (detailed balance), A |
| G8 | Generation | Impact ionization (vOdM), local Kane and nonlocal BTBT | physics, assemble | M, E (uniform field), A (breakdown) |
| G9 | Lattice heat | ∇·(κ∇T) + H = c ∂T/∂t, Wachutka fluxes, Peltier/Thomson | assemble | E (Kirchhoff, Joule), D (energy identity) |
| G10 | Energy balance | carrier temperatures, θn and θ²n SG fluxes | assemble | E (homogeneous field), P |
| G11 | Interfaces | heterojunction offsets, thermionic emission, Schottky (TE, image force), insulator interfaces, traps | assemble, device | E (equilibrium band diagram), A |
| G12 | Tunnelling | gate tunnelling (Tsu–Esaki, WKB), BTBT paths | physics, assemble | H (exact transmission), A (FN, Simmons), E (model integrals) |
| G13 | Quantum confinement | MLDA | physics | H (Schrödinger–Poisson), E (the MLDA factor) |
| G14 | Small signal | (J + iωC + T(iω)) δu = b | solve | E (linear devices), A (diffusion admittance), D (reciprocity, passivity) |
| G15 | Terminal quantities | weighted current, displacement current, charge | solve | D (KCL), E |

#### 2. Independence of the reference system (R2)

**2.1 Independence grades.** A reference is only as strong as it is independent of what it checks:

| Grade | Meaning | Accepted as a reference? |
|---|---|---|
| I3 | Different method (e.g. spectral collocation or shooting vs finite volume), separate code, different arithmetic (double-double), derived from the cited source | Yes: the target for E and H references |
| I2 | Different method, separate code, same arithmetic | Yes |
| I1 | The same formula re-coded separately from the source (no shared code) | Yes for M and A; for E only with a second I1/I2 cross-check |
| I0 | Calls NiTCAD's own physics or assembly code (a test that compares the code with itself) | **No.** Allowed as a regression pin, never as verification evidence |

**2.2 Independence rules (TVD-1, TVD-2 below).** The reference library depends on the C++ standard library and
`base`'s constants only (CODATA 2018, verified against NIST, ARCHITECTURE 10, V7); it never includes `physics`,
`assemble`, `solve`, `mesh`, `device` or `linalg` headers, enforced as `linalg`'s header-boundary test already is.
Material parameters reach both sides from the benchmark's own specification, never from NiTCAD's material tables, so
a wrong table entry is caught by L1, not hidden by sharing. Each reference function cites its source's equation
numbers. Each reference is reviewed as its own change, not with the code it checks.

**2.3 Audit of what exists.** A first task (TVR-0) classifies NiTCAD's existing analytic checks (about 270 matches for
analytic, exact, closed-form, Richardson and order terms across 55 test files, measured 2026-10-31) by class and
independence grade. For example, the MOS-C tests take their landmarks from the ported legacy model
(`tests/solve/legacy_moscap.hpp`): an independent but legacy-derived reference (I1 at best), useful, and to be
recorded as such. Nothing is removed or changed; the audit only labels and links them in the catalogue.

#### 3. Independent reference calculation architecture

**3.1 Layers (new; nothing existing changes):**

```text
verify/reference  : base only. Numerics (double-double, quadrature, roots, ODE/BVP, spectral collocation),
                    reference models (E, D, A, M, H), each with assumptions, applicability and its own uncertainty.
verify/consistency: results only. Conservation, equilibrium, reciprocity and passivity checks over any result.
verify/study      : drives the public solve API: mesh, time-step and tolerance ladders; Richardson and GCI.
verify/benchmark  : the catalogue (§6): device spec + analysis + reference binding + metrics + criteria + tier.
verify/report     : a versioned, checksummed verification report (io conventions, Unit 42 provenance) and CSV.
cli               : `nitcad verify [--tier fast|full|nightly] [--id TVB-…] [--report <file>]`.
```

**3.2 Reference numerics (each verified at L0 before use).** MSVC's `long double` is the 64-bit double, so higher
precision needs a **double-double** type (about 32 significant digits; Dekker/Knuth error-free transformations);
**adaptive Gauss–Kronrod (7/15)** quadrature with error estimate, and **tanh–sinh** for endpoint singularities (the
Poisson first-integral quadratures near the depletion edge); **bracketed Brent** root finding; **adaptive
Runge–Kutta (Dormand–Prince)** with dense output for shooting; **Chebyshev collocation** with domain splitting at
junctions and interfaces for 1D boundary-value problems (spectral convergence on each smooth piece: an I3 method
against the box method); **complex-step differentiation** to verify every manufactured source term and every
analytic derivative inside the reference. Each primitive has its own L0 tests on problems with known answers and
reports an error estimate, which becomes the reference's uncertainty U_ref.

**3.3 The reference model contract (R3, R10).** Every reference model declares, as data the report prints:
- *identity*: id, class (§1.2), independence grade (§2.1), citation with equation numbers, version;
- *equations*: the subset of G1–G15 and the exact form solved;
- *assumptions*: listed, each either built into the benchmark geometry (1D, uniform doping, abrupt junction) or a
  checkable condition;
- *applicability*: machine-checked predicates on the benchmark's parameters and on the solved state (e.g. low
  injection Δn/N ≤ 0.01 checked on the solution, not assumed; depletion approximation V_bi − V ≫ V_T; gradual channel
  L ≫ depletion depths; WKB κd ≫ 1). **A comparison outside its applicability is refused, not failed** (TVD-6);
- *outputs*: named quantities with units and U_ref (the reference's own numerical uncertainty);
- *known approximation error* (A class): an analytic bound or an estimate, with its source, or "measured" when an E
  reference for the same case exists and the difference E − A is reported;
- **provenance, derivation and verification (owner condition 3, TVD-18):** (i) *provenance*: the primary source
  (authors, journal, year, equation and table numbers) and how it was accessed, graded **PA** (primary read in full),
  **PB** (primary inaccessible; two or more independent reproductions read and agreeing), **PC** (one secondary
  reproduction), **PD** (derived here from the governing equations, no external source); (ii) *independent
  derivation*: the steps from the governing equations (G1–G15) to the formula used, written out in the reference's
  documentation, by someone other than the author of the code it checks where possible, with every approximation named;
  (iii) *verification*: L0 evidence (a second method agreeing, the source's printed values reproduced, limits checked).
  An **E** reference must be PA or PD-with-two-methods; an **A** or **M** reference may be PB with the grade printed in
  every report line; **PC is never enough on its own** (it is marked "provenance insufficient" and cannot pass a
  criterion).

**3.4 Extensibility (R10).** Adding a theory is: one reference model (contract above) in `verify/reference`, its L0
tests, and catalogue entries; nothing else changes. A checklist (TVD-12) gates each addition: source read and cited,
class and independence grade justified, applicability predicates tested at their boundary, U_ref measured, an
independent cross-check where class E. Future models already foreseen: density gradient (H: Schrödinger–Poisson), full
hydrodynamic (P), optical generation (E: Beer–Lambert slab), mixed-mode circuits (E: linear networks), AMR (Unit 17's
AD11 reuses §4's Richardson machinery).

#### 4. Comparison, error metrics and acceptance (R4, R5, R8)

**4.1 Error decomposition.** For a simulated quantity q_sim and reference q_ref:
q_sim − q_ref = e_disc + e_iter + e_lin + e_round + e_ref + e_approx + e_param, where e_disc is discretization
(mesh, time step), e_iter the Newton stopping error, e_lin the linear solve, e_round rounding (the current resolution,
Unit 15/32), e_ref the reference's own error (U_ref), e_approx the approximation error of an A or H reference, e_param
any parameter mismatch (zero by construction, §2.2). The framework estimates each (§4.3) and reports them beside the
discrepancy, so that a failure points at a term.

**The terms are not assumed additive or independent (owner condition 2, TVD-17).** The decomposition above is a
bookkeeping identity only where each term is defined as a difference between two computations that differ in one
factor; it is not a statement that the factors act independently. They interact: the Newton stopping error depends on
the mesh (the residual's scale changes with h); rounding grows as h shrinks (the current resolution of Units 15 and 32);
an A-reference's approximation error depends on the regime the mesh resolves; the observed order depends on whether
iterative error is below discretization error. Hence: (a) each term is estimated *at the state actually compared*, never
transplanted from another mesh or tolerance; (b) interactions are measured, not assumed away, by a factorial ladder on
the factors that matter for the benchmark (at least the 2×2 design h ∈ {h, h/2} × tol ∈ {tol, tol/100}; time step added
for transients), the interaction effect being e(h/2, tol/100) − e(h/2, tol) − e(h, tol/100) + e(h, tol), reported as
ε_int; (c) the combined bound is U = Σ_k \|U_k\| + \|ε_int\| (a conservative sum of magnitudes, not a root-sum-square,
which would assume independence), and a benchmark whose ε_int exceeds 10% of U reports "interacting errors" and cannot
pass an E criterion until the dominant interaction is resolved (e.g. the tolerance tightened until e_iter ≪ e_disc).

**4.2 Metrics (TVD-7):**

| Metric | Definition | Used for |
|---|---|---|
| Relative error with floor | \|q_sim − q_ref\| / max(\|q_ref\|, q_floor), q_floor = the current resolution for currents | scalars |
| Log error | \|log10(q_sim / q_ref)\|, both above their floors | quantities spanning decades (leakage, tunnelling, subthreshold current) |
| Field norms | discrete L∞; L2 weighted by dual volumes; H1 seminorm Σ_e (A_e / d_e) (Δe)² over dual faces; densities compared in quasi-Fermi potential or log density, never raw | fields (ψ, φ_n, φ_p, T) |
| Observed order | p̂ = ln(e_k / e_{k+1}) / ln r on ≥ 3 meshes of a family with ratio r; with no exact solution, from three solutions | L2, L3 |
| Asymptotic-range check | successive p̂ within 0.1 of each other, or the GCI ratio within 10% of r^p̂ | L3 (a p̂ outside the asymptotic range is not reported as an order) |
| GCI (Roache) | F_s \|ε₁₂\| / (r^p − 1), F_s = 1.25 with ≥ 3 meshes, 3.0 with 2 | U_disc |
| Iterative error | sensitivity × last Newton correction (Unit 32's adjoint form for currents) | U_iter |
| Conservation residuals | KCL: \|Σ_k I_k\| / max_k \|I_k\|; Gauss: \|∮ D·n − Q\| / max(\|Q\|, ε₀ E_scale A); energy: \|Σ I V − P_heat\| / Σ\|I V\| | §5 |
| Reciprocity, passivity | \|Y_ij − Y_ji\| / max\|Y\|; the smallest eigenvalue of (Y + Yᴴ)/2 ≥ −tol | small signal at equilibrium |
| Validation metric (X) | the comparison error against the data's stated uncertainty (ASME V&V 20 form: E = S − D, u_val from u_num, u_input, u_D) | L5 (Unit 33 catalogue) |

**4.3 Numerical uncertainty U_num.** The components: U_disc (GCI), U_time (the same on a time-step ladder), U_iter,
U_lin (the linear solver's backward error times the condition estimate where available), U_round (current resolution);
combined by §4.1's rule, U_num = Σ \|U_k\| + \|ε_int\|, with ε_int from the factorial ladder (not a plain sum of
independently estimated terms).
Each benchmark reports all five; a benchmark without its ladder reports U_num as "not estimated", and cannot pass an
E or A criterion (TVD-8).

**4.4 Acceptance criteria by class (TVD-9; fixed per benchmark in the catalogue, never loosened after a run):**

| Class | Criterion |
|---|---|
| D | Holds to rounding: residual ≤ c·ε·(condition-scaled magnitude), c stated per benchmark (default 64, Unit 50's rounding-floor convention). |
| E | (a) **Convergence:** p̂ within [p_th − 0.2, p_th + 0.3] in the asymptotic range on ≥ 3 meshes, where p_th is the order the scheme is proven or established to have for that norm (§5.3); (b) **accuracy:** \|q_sim − q_ref\| ≤ U_num + U_ref on every mesh and ≤ the benchmark's tolerance on the finest; (c) U_ref ≤ tolerance / 10. |
| A | Inside the applicability predicates only: \|q_sim − q_ref\| ≤ U_num + U_ref + ε_approx, ε_approx from the reference's stated bound or from a measured E − A difference on the same case. A pass says "consistent with the theory within its approximation", never "exact". |
| M | The implementation equals the cited formula to 1e-13 relative over a parameter grid covering the model's stated range, and the source's tables to their printed digits; outside the stated range: reported, not accepted. |
| H | No pass/fail: the model-form error is measured and recorded against the regime (e.g. MLDA vs Schrödinger–Poisson inversion charge vs oxide field). A threshold may be set by an owner ruling. |
| P | Within the published plot's reading uncertainty; a miss is a finding. |
| X | Unit 33's existing pins, unchanged. |

**Outcomes (owner condition 1, TVD-16):** pass; **finding** (a measured, pinned discrepancy with its attribution and an
owner ruling, as Units 33, 48 and 54 did); fail; **NOT APPLICABLE** (an applicability predicate of the reference or
check is false for this device or state; the report names the predicate and its measured value, e.g. "C7 reciprocity:
NOT APPLICABLE, bias ≠ 0"); **NOT EVALUATED** (the inputs to decide applicability are missing, e.g. no ladder). NOT
APPLICABLE is never counted as a pass, and a benchmark whose every check is NOT APPLICABLE is reported as unverified.
**Rule TVM2 (priority over passing tests):** a criterion is never widened to make a run pass; a reference is never
adjusted toward the simulation; a model parameter is never tuned to a benchmark (Unit 33's V-rules carried over).

#### 5. Physical consistency, conservation and numerical convergence (R6)

**5.1 Consistency checks (`verify/consistency`, run on any result, not only benchmarks).** Each check carries its
applicability predicates (the "Holds when" column, evaluated on the device, models and state); where one is false the
outcome is NOT APPLICABLE with the predicate named (TVD-16), e.g. C1/C2 off equilibrium, C7 under bias or a magnetic
field, C8 with active generation, C9 outside the linear regime, C10 when the limit parameter is not swept.

| # | Check | Class | Holds when |
|---|---|---|---|
| C1 | Equilibrium: zero terminal and edge currents, flat quasi-Fermi levels, n p = n_ie² at every node | D | Zero bias, no generation (SG is exact here) |
| C2 | Detailed balance: every G − R term is zero at equilibrium (SRH, Auger, radiative, traps, II, BTBT) | D | Equilibrium |
| C3 | Discrete Gauss law: flux of D through every dual cell equals its charge; global charge neutrality including electrode charges | D | Any converged state, to the Newton tolerance |
| C4 | KCL: Σ terminal currents = 0 (DC); with displacement current (transient, AC) | D | Any state; the weighted current's resolution |
| C5 | Energy: Σ I V = ∫ H dV + boundary heat flow (electrothermal), carrier energy balance closure (Unit 52) | D | Steady states |
| C6 | Transient charge balance: ∫ I dt = ΔQ per terminal | D | Transient runs, per step and cumulative |
| C7 | Reciprocity (Onsager): Y(ω) symmetric at equilibrium for a device without magnetic field | D (exact for the linearized symmetric operator) / measured otherwise | AC at zero bias |
| C8 | Passivity: Re(Y) positive semidefinite at equilibrium | E (thermodynamics) | AC at zero bias |
| C9 | Time–frequency consistency: the AC admittance at low ω equals the DC conductance and the quasi-static capacitance; the transient small-step response's Fourier transform matches Y(ω) in band | E | Linear regime |
| C10 | Limit continuity: FD → Boltzmann (non-degenerate), energy balance → drift-diffusion (τ_w → 0), MLDA → classical (far from interfaces), Lombardi → bulk (E⊥ → 0), electrothermal → isothermal (κ → ∞) | E (limits of the models) | Parameter sweeps toward the limit; the difference must vanish at the stated rate |
| C11 | Symmetry and invariance: mirror-symmetric devices give mirrored fields; rigid rotation or translation of an unstructured mesh leaves terminal quantities unchanged to rounding; dimensional reduction (ARCHITECTURE 10's existing gate) | D | — |
| C12 | Positivity and bounds: n, p > 0; T > 0; occupancies in [0, 1] | D | Always |

Several exist inside unit tests today (C1, C4, C5 partly, C11's reduction); §5 makes them a reusable layer reported for
every benchmark.

**5.2 Numerical convergence studies (`verify/study`).** Mesh ladders (uniform refinement ratio r = 2 per family:
tensor, 2D generated, swept, 3D generated); time-step ladders for transients (expected orders 1 for backward Euler, 2
for BDF2); Newton-tolerance ladders (the iterative error must be ≤ 1/10 of U_disc before an order is computed);
linear-solver checks (backward error, and SparseLU vs PARDISO vs GMRES agreement to their stated accuracy).

**5.3 Expected orders (stated with their evidence, not assumed).** For the box method with Scharfetter–Gummel fluxes
on admissible (Delaunay–Voronoi) meshes: first order in the discrete H1 seminorm is stated to be proven for elliptic
problems (Eymard, Gallouët and Herbin's finite-volume analysis); second order in L2 and at nodes is typical on smooth
solutions and on uniform tensor meshes, but not guaranteed on general unstructured meshes. **[unverified, TVR-0
§10.6]:** the theorem's text, number and hypotheses have not been read (the open copy, hal-02100732, refused automated
access on 2026-10-31); these orders are working assumptions until it is read. Abrupt doping and material
interfaces resolved only to h reduce orders (Unit 52 measured first order at abrupt junctions). Hence: E criteria use
p_th = 1 for H1 on unstructured families, p_th = 2 for L2/L∞ on tensor families and on smooth manufactured solutions;
L2 orders on unstructured families are reported, and gated only where a family's order is first established by a
smooth benchmark (TVD-10).

#### 6. Validation benchmark catalogue (R4, R8)

Columns: id; what; equations (§1.3); class (§1.2); reference method and source; quantities and metric; criterion
(§4.4) and key assumptions/applicability; priority (P1 first). Tolerances marked "set at TVR-1" are fixed from the
reference's U_ref and the target mesh before any simulation is run, then frozen.

| ID | Benchmark | Eqs | Class | Reference (method; source) | Quantities, metric | Criterion; assumptions | Pri |
|---|---|---|---|---|---|---|---|
| TVB-01 | Constants and thermal voltage | — | E | CODATA 2018 | values, exact | D-type: bit-equal to the stated decimal | P1 |
| TVB-02 | Fermi–Dirac integrals F_{1/2}, F_{−1/2} and inverse | G4 | E | Double-double adaptive quadrature; Blakemore and Bednarczyk/Aymerich-Humet as cross-checks | relative error over η ∈ [−40, 40] | ≤ 1e-12 relative; defines L0 for statistics | P1 |
| TVB-03 | Bernoulli function B(x) = x/(eˣ − 1) near 0 and at large \|x\| | G3 | E | Double-double series and direct | relative error | ≤ 4ε | P1 |
| TVB-04 | Linear Poisson, manufactured (1D/2D/3D; tensor, 2D generated, 3D generated, swept) | G1 | E (MMS) | Smooth ψ with ε jump on a planar interface; sources by complex-step-checked derivatives | L2, H1, L∞, p̂ | §4.4 E with p_th per §5.3 | P1 |
| TVB-05 | Coaxial and spherical capacitors on unstructured meshes | G1 | E | C = 2πεL / ln(b/a); C = 4πε ab/(b − a) | capacitance, p̂ | E; faceting error O(h²) included in p_th | P1 |
| TVB-06 | Dielectric stack (MOS insulator, high-k) | G1 | E | Series capacitance, exact | C, field ratio ε₁/ε₂ | D-type on aligned meshes | P1 |
| TVB-07 | Equilibrium abrupt pn junction, Boltzmann | G1, G4 | E (semi-analytic) | First integral of 1D Poisson–Boltzmann: junction field in closed form, x(ψ) by tanh–sinh quadrature, ψ_j by Brent; I3 cross-check by Chebyshev collocation | V_bi (exact: V_T ln(N_A N_D / n_i²)), E_max, ψ(x), depletion edges; L∞, p̂ | E; uniform doping each side, constant ε | P1 |
| TVB-08 | Equilibrium heterojunction (Anderson rule), with ε step | G1, G11 | E (semi-analytic) | First integrals on each side matched by D continuity | band diagram, spike/notch heights, E on each side | E; abrupt, no interface charge | P2 |
| TVB-09 | MOS capacitor, uniform substrate: Q_s(ψ_s), LF C-V | G1, G4 | E (semi-analytic) | Kingston–Neustadter F-function (exact 1D); V_G = V_FB + ψ_s − Q_s / C_ox; C_LF = dQ/dV_G by complex step | Q_s, C_LF(V_G), V_FB, V_T landmark; relative error, p̂ | E; Boltzmann, uniform doping, 1D | P1 |
| TVB-10 | MOS-C HF C-V and deep depletion | G1, G14 | A | Minority charge frozen (HF approximation) on the exact 1D solution | C_HF, C_min | A; HF approximation (ω ≫ generation rate) | P2 |
| TVB-11 | Uniform resistor (ohmic contacts) | G1–G3 | D/E | V = I L / (q μ N A); SG exact at constant field | I(V), ψ linear; to rounding | D | P1 |
| TVB-12 | Scharfetter–Gummel constant-field exactness | G3 | D | 1D, constant E and J, no G − R: nodal densities equal the exact exponential solution | nodal n; to rounding | D | P1 |
| TVB-13 | Uniform-field drift with Canali velocity saturation; homogeneous energy balance | G6, G10 | E | v(E) = μ₀E / (1 + (μ₀E/v_s)^β)^{1/β}; energy balance: algebraic T_n(E) of the homogeneous steady state | v(E), T_n(E); relative | E (model's own homogeneous limit) | P2 |
| TVB-14 | Full drift-diffusion, manufactured (ψ, n, p; 1D/2D/3D; Boltzmann, SRH on) | G1–G3, G7 | E (MMS) | Manufactured positive densities (exponential forms); sources by complex step in double-double; **needs a source-term hook** (TVD-15) | L2/H1 in ψ, φ_n, φ_p; p̂ | E | P1 |
| TVB-15 | Neutral-region minority diffusion with SRH (low injection) | G2 | E (linearised) + A | Exact solution of the linear minority ODE (cosh/sinh); low injection checked on the solved state | minority profile, injected current | E for the linear model; A vs full DD with ε_approx = measured Δn/N | P2 |
| TVB-16 | Ideal diode: long and short base; Sah–Noyce–Shockley recombination current | G1–G3, G7 | A | Shockley (J₀ = q n_i² Σ D/(L N)); SNS integral | J(V), ideality; log error | A; low injection, V_bi − V ≫ V_T; the legacy 1.321e-2 short-base value as cross-check | P1 |
| TVB-17 | Depletion capacitance | G1, G14 | E (semi-analytic) + A | Exact dQ/dV from TVB-07's first integral at reverse bias; ε/W with the −2V_T correction as the A form | C_j(V) | E vs exact; A difference measured | P2 |
| TVB-18 | Long-diode small-signal admittance | G2, G14 | A | Y = G₀ √(1 + iωτ) + iωC_j | Y(ω) magnitude and phase | A; low injection, long base | P2 |
| TVB-19 | Dielectric relaxation and Haynes–Shockley pulse | G2 (transient) | E (linear) | τ_d = ε/σ; drifting decaying Gaussian for the linear minority equation in a uniform field; **needs an initial excess-density state via the transient API** (exists: initial state) | decay rate, pulse centre, width; time-step p̂ | E with BDF2 p_th = 2 | P2 |
| TVB-20 | Long-channel MOSFET: Pao–Sah double integral; charge sheet; threshold and body effect | G1–G3, G6 | A (Pao–Sah exact under gradual channel) | Pao–Sah by nested adaptive quadrature over TVB-09's exact 1D charge; Brews charge sheet as second A | I_D(V_G, V_D), V_T, g_m | A; constant mobility, L ≫ depletion depths (checked) | P2 |
| TVB-21 | Subthreshold swing | G1–G3 | A | S = ln10 V_T (1 + C_dep/C_ox) (+ C_it) | S | A; long channel, low V_D | P2 |
| TVB-22 | Interface traps: single-level conductance and capacitance | G7, G11, G14 | E (model) | Nicollian–Brews equivalent circuit for one level (exact for the small-signal SRH trap model) | G_p/ω peak, C_it | E | P3 |
| TVB-23 | Impact ionization: avalanche breakdown and multiplication | G8 | A (+ E for frozen field) | Ionization integral with vOdM coefficients over the exact 1D field; depletion-field version as A | V_BR, M(V) | E for the frozen-field integral; A for V_BR | P3 |
| TVB-24 | BTBT: Kane rate in a uniform field; nonlocal path integral | G8 | E (model) | Closed-form Kane rate; path integral by quadrature | G, current | M and E | P3 |
| TVB-25 | Gate tunnelling: trapezoidal barrier | G12 | H + A | Exact transmission with Airy functions (Gundlach) inside the Tsu–Esaki integral; WKB and Fowler–Nordheim and Simmons as A | J(V_ox); log error | H: WKB model-form error measured by regime; A: FN in its regime | P3 |
| TVB-26 | Schottky diode: thermionic emission, image-force lowering; thermionic-diffusion | G11 | A | J = A** T² e^{−φ_B/V_T}(e^{V/V_T} − 1); Crowell–Sze combination | J(V), ideality, barrier lowering | A; Unit 48's accepted limitations recorded beside | P3 |
| TVB-27 | Lattice heat: Kirchhoff slab (T-dependent κ); Joule-heated resistor | G9 | E | Kirchhoff transform (existing U4 gate, re-homed); parabolic T for constant κ | T(x), p̂ | E | P2 |
| TVB-28 | Thermoelectric identities: Seebeck, Peltier, energy identity | G9 | E (thermodynamic) | Existing DC7 values re-homed | S, Π, energy residual | E/D | P2 |
| TVB-29 | Cylindrical and spherical pn junctions (curved 3D unstructured) | G1, G4 | E (semi-analytic) | Radial Poisson–Boltzmann by Chebyshev collocation | V_bi, E_max, C_j | E | P3 |
| TVB-30 | MLDA vs 1D Schrödinger–Poisson (effective mass, infinite and finite barrier) | G13 | H | Independent 1D Schrödinger–Poisson: the finite-difference Hamiltonian is symmetric tridiagonal, solved by an eigensolver written in the reference (bisection on Sturm sequences for the eigenvalues, inverse iteration for the vectors, in double-double; I3), with Eigen's `SelfAdjointEigenSolver` as a secondary cross-check only (owner ruling on TVQ4) | inversion charge centroid, ΔV_T, C_g | H: model-form error recorded vs field and doping | P3 |
| TVB-31 | Model-implementation catalogue: Caughey–Thomas, Masetti, Canali, Lombardi, Slotboom, Couderc band set, Auger, radiative, vOdM, Kane, MLDA factor, Tsu–Esaki supply function | G4–G8, G12, G13 | M | Formulas re-coded from the cited papers (I1); source tables | relative error over parameter grids | M (1e-13; tables to printed digits) | P1 |
| TVB-32 | Consistency suite C1–C12 on every P1/P2 benchmark device | all | D/E | §5.1 | residuals | D/E per check | P1 |
| TVB-33 | Published numerical: Aste–Vahldieck ballistic diode (energy balance) | G10 | P | The paper's plotted peak velocity | peak v | P (existing Unit 52 measurement re-homed) | P3 |
| TVB-34 | Published data: Unit 33's resistivity and n_i benchmarks | — | X | Thurber et al.; Misiakos–Tsamakis | Unit 33's pins | X, unchanged | P1 (catalogue only) |

Not in the catalogue now, with the reason: velocity overshoot and full hydrodynamics (no accepted closed form; P class
only); FinFET-scale devices (no reference beyond tensor self-convergence: Unit 54's UC8 is L3, not L4); breakdown
from measured data (no accessible source, Units 36–39 plan); optical generation, circuits (NiTCAD lacks them).

#### 7. Limitations and sources of discrepancy (R9)

- **The reference can be wrong.** Mitigations: L0 tests for every numerical primitive; class E references
  cross-checked by a second independent method (TVB-07, TVB-09, TVB-29); hand values from the source's own tables
  wherever printed; separate review (§2.2). A disagreement is investigated on both sides before it is attributed.
- **Approximation theories are not truth.** A-class comparisons measure consistency within an approximation;
  ε_approx is itself an estimate unless an E reference for the same case exists.
- **Empirical models.** M-class checks verify transcription, never physics; agreement says nothing about accuracy
  against devices (that is L5).
- **Discrete vs continuous models.** Abrupt features resolved only to h, interfaces represented by straddles,
  faceted curved surfaces, and lumped oxides all change the continuous model the discretization approximates; the
  catalogue states for each benchmark which continuous model is meant.
- **Rounding floors.** Minority densities near 1e-35 (4H-SiC), currents below the resolution (Units 15, 32, 50): log
  metrics and floors are mandatory; claims below the floor are refused.
- **Nonlinear multiplicity.** Near folds (breakdown, snap-back, thermal runaway) the solution branch matters;
  benchmarks there compare on the traced branch (Unit 19's arc length) and name it.
- **Platform effects.** Last-bit differences across P/E cores were found (Unit 40, fixed by Unit 41) and across
  backends; D-class tolerances are stated relative to rounding, not bit equality, unless the benchmark pins a backend.
- **Hooks NiTCAD lacks.** Manufactured-solution sources (TVB-14) need an equation-term hook (Unit 40's generic terms
  make it a small, test-only term); prescribed generation for optical cases is absent.
- **Cost.** Mesh ladders on 3D are expensive (Unit 53's FC12); 3D ladders go to the nightly tier.
- **Coverage is never complete.** The catalogue covers the equations in §1.3, not every combination of models; a
  combination without a benchmark is reported as unverified, not as verified by its parts.

#### 8. Design choices and owner questions

**Design (frozen for review; no code):**

| # | Choice | Why |
|---|---|---|
| TVD-1 | `verify/reference` depends on the standard library and `base` constants only; a header-boundary test enforces it. | Independence (§2) |
| TVD-2 | Parameters flow from the benchmark spec to both sides; never from NiTCAD's tables to the reference. | No shared-table masking |
| TVD-3 | Double-double, Gauss–Kronrod, tanh–sinh, Brent, Dormand–Prince, Chebyshev collocation, complex step, all written in the reference with their own L0 tests. | I3 references; no MSVC extended precision |
| TVD-4 | Every reference model implements the contract of §3.3 (identity, equations, assumptions, applicability, outputs with U_ref, approximation error). | R3, R7, R10 |
| TVD-5 | Benchmarks are data plus a device builder (the deck where expressible; the API where a hook is needed); each states tier, class, criterion and its frozen tolerance. | Reproducible, reviewable |
| TVD-6 | Applicability predicates are evaluated on the parameters and the solved state; out-of-range comparisons are refused with the predicate named. | Rigour over passing |
| TVD-7 | Metrics of §4.2, in one library; densities compared in log or quasi-Fermi form only. | Comparable across benchmarks |
| TVD-8 | No pass for E or A without U_num (§4.3). | An unconverged agreement is luck |
| TVD-9 | Criteria per class (§4.4), frozen before the first run; outcomes pass, finding, fail. | Owner's priority on rigour |
| TVD-10 | Expected orders per family from §5.3, with the evidence recorded. | Orders are claims, not hopes |
| TVD-11 | The verification report: versioned, checksummed, with NiTCAD's provenance (Unit 42), each reference's version, U_ref, U_num breakdown, the outcome and attribution; two runs bit-identical on one machine and thread count. | Reproducibility (R8) |
| TVD-12 | An addition checklist for new references (§3.4). | Extensibility (R10) |
| TVD-13 | Tiers: `fast` (D, M, small E; Debug-safe, ≤ 60 s added), `full` (Release, all P1–P2), `nightly` (3D ladders, P3). Release-only like Unit 33's `validation` ctest. | CI cost |
| TVD-14 | Existing tests are not moved or changed; the audit (TVR-0) labels them and links them from the catalogue; re-homing a check into `verify` happens only with the original gate kept. | "Do not modify existing architecture" |
| TVD-15 | Dependencies listed, not designed here: a test-only manufactured-source equation term (TVB-14), transient initial excess state (exists; TVB-19). | Scope |
| TVD-16 | Outcomes include NOT APPLICABLE (predicate named, measured value printed) and NOT EVALUATED; neither counts as a pass (owner condition 1). | Assumptions are checked, not presumed |
| TVD-17 | Error terms estimated at the compared state; interactions measured by a factorial ladder; combined as a sum of magnitudes plus the measured interaction; large interactions block E passes (owner condition 2). | No assumed additivity |
| TVD-18 | Every reference carries provenance grade PA–PD, a written independent derivation and L0 verification; E references PA or PD with two methods; PC insufficient (owner condition 3). | Traceable references |

**Owner questions (ruled 2026-10-31: all approved with conditions 1–3; TVQ4 amended to an in-reference symmetric
tridiagonal eigensolver with Eigen as a secondary cross-check; see the rulings at the top of this entry):**

| # | Question | Recommendation |
|---|---|---|
| TVQ1 | Accept the class scheme (E, D, A, M, H, P, X) and independence grades (I0 not accepted as evidence)? | Yes |
| TVQ2 | Accept that existing analytic checks are audited and labelled, not rewritten (TVD-14)? | Yes |
| TVQ3 | The manufactured-solution hook (TVB-14): a test-only equation term in `assemble`, behind the generic terms of Unit 40, bit-identical when absent. It is the only change to existing code the framework needs for its P1 set. | Approve as part of phase 1 |
| TVQ4 | Use Eigen in the reference for dense eigenproblems (TVB-30, grade I2), or write a reference eigensolver (I3)? | Eigen (third-party, independent of NiTCAD's physics); revisit if TVB-30 becomes a gate |
| TVQ5 | H-class (model-form) comparisons: record only, or gate with owner-set thresholds? | Record only in phase 3; thresholds by later ruling |
| TVQ6 | The roadmap's split into implementation units (below)? | As proposed |

#### 9. Prioritised implementation roadmap (each phase its own unit, proposal and criteria before code)

| Phase | Contents | Benchmarks | Size |
|---|---|---|---|
| TVR-0 | Audit of existing analytic checks: class, grade, catalogue links (documentation only) | — | S |
| TVR-1 (Unit 55a) | `verify/reference` numerics with L0 tests; metrics; contract; report format; `nitcad verify`; consistency suite C1–C12 | TVB-01–03, 32 | L |
| TVR-2 (Unit 55b) | Code verification core: MMS Poisson and drift-diffusion (with the TVQ3 hook), SG exactness, resistor, capacitors, equilibrium pn junction, MOS-C exact; convergence study driver with GCI | TVB-04–07, 09, 11, 12, 14 | L–XL |
| TVR-3 (Unit 55c) | Model-implementation catalogue (M class) for every model in NiTCAD | TVB-31 | M |
| TVR-4 (Unit 55d) | Device theory (A class) and small signal / transient: diode family, depletion capacitance, admittance, Pao–Sah, swing, relaxation, Haynes–Shockley, thermal | TVB-10, 13, 15–21, 27, 28 | L |
| TVR-5 (Unit 55e) | Advanced physics: heterojunction, traps, II, BTBT, tunnelling (exact Airy), Schottky, curved junctions, MLDA vs Schrödinger–Poisson, published numerical | TVB-08, 22–26, 29, 30, 33 | XL |
| TVR-6 | Integration: tiers in CI, the report in the application, Unit 17's AD11 reusing §4's machinery, Unit 33's catalogue linked | TVB-34 | M |

**Order rationale:** L0 (the reference's own correctness) before anything it judges; code verification (L2) before
theory comparison (L4), per TVM1; M-class early because it is cheap and catches parameter transcription before any
device-level attribution. Unit 17 (adaptive refinement) benefits from TVR-1/TVR-2 (its certificate's Richardson
estimate and its references); the two can proceed independently, and TVR-2's benchmarks become Unit 17's accuracy
references when both exist.

**Sources named (to be read and cited by equation number at each phase's proposal; not all read yet):** Roache,
*Verification and Validation in Computational Science and Engineering* (1998); Oberkampf and Roy, *Verification and
Validation in Scientific Computing* (2010); ASME V&V 20-2009; Salari and Knupp, *Code Verification by the Method of
Manufactured Solutions* (SAND2000-1444); Eymard, Gallouët and Herbin, *Finite Volume Methods* (Handbook of Numerical
Analysis VII, 2000); Sze and Ng, *Physics of Semiconductor Devices* (3rd ed.); Pierret, *Advanced Semiconductor
Fundamentals*; Nicollian and Brews, *MOS Physics and Technology*; Pao and Sah (Solid-State Electronics 9, 1966);
Brews (Solid-State Electronics 21, 1978); Gundlach (Solid-State Electronics 9, 1966); Simmons (J. Appl. Phys. 34,
1963); Crowell and Sze (Solid-State Electronics 9, 1966); Blakemore (Solid-State Electronics 25, 1982); Dekker
(Numer. Math. 18, 1971); and the sources of each M-class model as already cited in DECISIONS.md.

#### 10. TVR-0: audit of the existing validation (research; 2026-10-31 to 2026-11-04; nothing changed, nothing committed)

**10.1 Scope and method.** The suite has 1,163 test cases; the 627 in `tests/physics`, `assemble`, `solve`, `analysis`,
`validation` and `device` were listed (name and line), and every one whose name or header claims a theoretical,
analytic, exact, published or convergence reference was read at its reference computation (about 90 cases, in
`equilibrium`, `bias`, `heterojunction`, `insulator`, `small_signal`, `impact_ionization`, `electrothermal`,
`energy_balance`, `unstructured`, `transient`, `weighted_current`, `quantum_confinement`, `mosfet`, `deck/schottky` and
the `physics` model tests), together with the shared reference helpers `tests/physics/references.hpp` and
`reference_arithmetic.hpp` and the legacy ports `tests/solve/legacy_*.hpp`. The rest (Jacobian, refusal, record,
cost, cancellation and I/O cases) were classified from their names and headers only; they are discrete-consistency or
software checks, not theory comparisons. `app`, `render`, `io`, `deck` vocabulary, `linalg` and `mesh` were out of
scope except the mesh convergence gates. Independence grades are those of §2.1; classes those of §1.2; provenance grades
those of §3.3.

**10.2 What exists, by kind.**

| Kind | Representative cases (file:line) | Class / grade | Assessment |
|---|---|---|---|
| K1 Discrete consistency | FD-Jacobian gates (every `assemble` file); zero current at equilibrium (`bias:318`, `heterojunction:191`); KCL and Gauss (`mos:343`, `insulator:382`, `transient:365`, `transient:533`); SG constant current (`bias:331`); Laplace exact on a Voronoi dual (`unstructured:69`); y-uniform reduction to 1D (about 20 cases); reciprocity at equilibrium (`small_signal:341`); Y(0) = DC slope (`small_signal:369`); energy identities (`electrothermal:547`, `energy_balance:712`); div J = qU (`rates:171`) | D; I0 by nature | **Sound and strong.** Correctly framed as consistency, tolerances tied to rounding or the current resolution. Missing: passivity (C8), Kramers–Kronig (C9), limit continuity as a systematic check (C10, present only for EB → DD, `energy_balance:327`). |
| K2 Exact with an independent formula | Exact MOS-C first integral with traps and fixed charge, second order shown (`insulator:256`, `:287`, `:325`); heterojunction first integral (`heterojunction:147`); Kirchhoff transform (`electrothermal:138`, `:1013`); step-heated rod series (`electrothermal:1334`); self-heated resistor ODE (`electrothermal:650`); Seebeck/Peltier integrals (`:387`, `:417`, `:743`, `:783`); uniform-field EB, BTBT, II (`assemble/energy_balance:812`, `assemble/band_to_band:285`, `assemble/impact_ionization:214`); Fermi–Dirac integrals and Bernoulli in double-double (`physics/fermi_dirac:69`, `assemble/flux:27`) | E; formula I1–I3, **inputs often I0** | **The strongest part of the suite.** The double-double library (`reference_arithmetic.hpp`) was checked here: two_sum/two_prod by FMA, exp by ln 2 and 2⁻¹⁰ reduction with a 12-term Taylor series and ten squarings, log by Newton on exp, Gauss–Legendre by Newton on Pₙ, Bernoulli numbers B₂–B₂₀ and the Euler–Maclaurin ζ: all correct as written; the stated ~1e-30 accuracy is consistent with the 2¹⁰ squaring amplification of a ~1e-32 series. Weakness: the "exact" side often takes n_i, N_c, mobilities or ionization coefficients from NiTCAD's own physics functions (15 of the `solve` test files do), so an error in those functions cancels (C-9 below). |
| K3 Approximate theory | Depletion-approximation peak field (`equilibrium:156`); ideal-diode law (`bias:165`); reverse-leakage exponent (`bias:284`); depletion capacitance and C–V doping profile (`small_signal:856`, `:878`); Shockley diffusion admittance (`small_signal:994`); single-level conductance method (`small_signal:735`); one-sided breakdown (`impact_ionization:203`); long-channel swing bracket (`mosfet:222`); surface-recombination diode (`insulator:602`); trap emission rate (`transient:385`); Schottky SC3, SC5 (`deck/schottky:198`, `:362`) | A; I1 | Mostly correctly reasoned, but approximation errors are asserted rather than measured, tolerances are often far wider than the numerical error, and two are mis-classified (I-1, I-5 below). |
| K4 Empirical models | Masetti (`physics/masetti`), Couderc band set (`band_set`), Altermatt levels (`ionization_density`), Lombardi (`surface_mobility`), van Overstraeten–de Man (`physics/impact_ionization`), Canali, Caughey–Thomas, Slotboom, Kane, MLDA factor; Thurber fits (`validation/published:160`) | M; I1 | Formulas re-coded from DECISIONS.md, not from the implementation's structs (good). **Provenance is often secondary** (C-10). |
| K5 Self-convergence | Unstructured 2D/3D orders against a 4× finer tensor solve (`unstructured:119`, `unstructured3d:302`); abrupt diode against NiTCAD's 1D solve (`unstructured:200`); MOSFET rule against the tensor MOSFET (`unstructured:356`); FinFET/complex convergence (`complex_convergence:396`, `:452`); time orders against NiTCAD's own fine run (`transient:209`, `electrothermal:1334` time part); `bias:257`; `mosfet:499`; QC8; EC9 | L3 solution verification; I0 | Legitimate as solution verification, but **not code verification**: a defect common to both discretizations (same SG fluxes, same physics) is invisible (C-5). No GCI safety factor or asymptotic-range check anywhere (I-7). |
| K6 Legacy / cross-code regression | J(0.5 V) = 1.280e-2 (`bias:144`); legacy MOS-C port and P1–P9 (`mos:172`, `mos:328`–`:584`, `legacy_moscap.hpp`); legacy turn-off (`transient:718`); legacy G-PARABOLA (`electrothermal:602`); legacy MOSFET bands (`mosfet:222`) | Regression; "I0-lineage" | Useful continuity pins; **not theoretical evidence** (C-7, C-8). |
| K7 Published numerical / measured | Aste–Vahldieck ballistic diode (`energy_balance:375`); Unit 33 benchmarks (`validation/published`) | P; X | EC10 compares different models (I-8). Unit 33 is a careful X baseline (out of TVR-0's scope beyond provenance). |

**No manufactured-solution test exists** (no file mentions one): there is no code verification of the coupled
drift-diffusion system against an exact solution on unstructured meshes (§10.5).

**10.3 Self-referencing or circular validations.** Severity: **H** (a check presented as verification that cannot
detect the error class it names), **M** (partly circular, mitigated elsewhere), **L** (labelling only).

| # | Case | Circularity | Sev. | Mitigation present? |
|---|---|---|---|---|
| C-1 | `equilibrium:92` "built-in potential within 2e-3 V (legacy gate)", listed in ARCHITECTURE 10 as a validation criterion | Measures ψ(last) − ψ(first), two Dirichlet contact values: it tests the contact formula, not the solve (the file's own header says so, and so does the ARCHITECTURE 10 row); and the expected V_bi is built from NiTCAD's `effective_intrinsic_density` (`built_in()`, `equilibrium:83`), so a wrong n_ie passes | M (the n_ie circularity; the not-a-solve-gate limitation is already disclosed) | The solve is checked by `equilibrium:130` and `:156`; the ARCH 10 row already states that it "does not test the solve" and points to the peak-field gate |
| C-2 | `impact_ionization:203` "the primary check": breakdown where the ionization integral on the *simulated* field reaches 1; `impact_ionization:170` current balance | Both integrate NiTCAD's own α(E) over NiTCAD's own field and densities: they verify that the solver's runaway is the model's own criterion (a discrete consistency), not that the field or α is right | M | The analytic one-sided V_BR (A class) is the independent part, also with NiTCAD's α; α itself is checked against the paper's fit magnitudes (`physics/impact_ionization:22`) |
| C-3 | Schottky SC4 "series law on the solution's own potential" (`deck/schottky:286`) | The reference is evaluated on the computed potential | M | SC3 and SC5 compare with closed forms |
| C-4 | `small_signal:232`, `:271` (RC), `:681` (HF C), `:705` (gated diode LF C), `:783` (transient vs AC) | The "expected" C is NiTCAD's own quasi-static dQ/dV, or `frozen_minority_capacitance` re-solved on NiTCAD's DC state; R uses NiTCAD's `masetti_mobility` | M | These are valuable *inter-analysis consistency* checks (C9-type) and should be labelled so; none is a theory comparison for C |
| C-5 | All K5 self-convergence gates | The reference is NiTCAD's own finer (tensor or 1D) solve with the same fluxes and physics | M | `insulator:256` (exact MOS-C) and `heterojunction:147` are true exact comparisons, but in 1D only |
| C-6 | `transient:209`, `electrothermal:1334` time orders | Against NiTCAD's own 5120- or 1280-step run | L | `electrothermal:1334`'s spatial part is against the exact series |
| C-7 | `bias:144` J(0.5 V) = 1.280e-2 within 1% (ARCHITECTURE 10) | The number is the legacy pytcad simulator's output for the same algorithm: a cross-version regression pin. The test notes the short-base analytic value 1.321e-2 (3.2% away) but does not compare with it | M | `bias:165` compares with the analytic law (I-2) |
| C-8 | `mos:172` and the P1–P9 MOS-C gates | A port of the legacy MOS-C (independent residual and Newton, but the legacy's algorithm, mesh and shared material functions); P1 checks the legacy's own discrete rows | L | `insulator:256` (exact) covers the same physics independently |
| C-9 | References whose inputs come from the code under test: `physics::intrinsic_density`, `effective_intrinsic_density`, `caughey_thomas_mobility`, `masetti_mobility`, `impact_ionization_coefficient`, `conduction_band_dos`, `band_gap_eV` inside expected values (15 `solve` files; most in `heterojunction`, `insulator`, `impact_ionization`, `electrothermal`, `small_signal`); and `references.hpp`'s material formulas, evaluated in double-double "from this code's parameter sets" | An error in a parameter table or a derived-quantity function appears on both sides and cancels | M | L1 tests of those functions exist (`semiconductor:64` n_i, `band_set:72`, `mobility:18`, `masetti`), but several take their parameters from the code's own structs too (`mobility_test`, `heavy_doping_test` via `references.hpp`), so a wrong *parameter value* in `legacy_silicon` is caught only where a literal from a source is compared (`mobility:18`, `semiconductor:64`) |
| C-10 | M-class provenance (§3.3 grades): Masetti 1983 (paywalled; three reproductions, **PB**); Masetti phosphorus μ_max = 1414 (one reproduction, **PC**, DECISIONS Unit 38 limitations); Altermatt 2006 (paywalled, not read, **PB/PC**); Lombardi 1988 (paywalled; two text reproductions, equations as images not read, **PB**); MLDA (Paasch–Übensee, paywalled; density equation an image, the Boltzmann reduction derived here, **PD** for the reduction); Schottky references (equations as images not read); Thurber (OCR, some coefficients reconstructed, DECISIONS Unit 33); Sze–Gibbons breakdown (not read, V6) | The tests verify agreement with DECISIONS.md's transcription of secondary reproductions, not with the primary source | M (PC: H) | Recorded honestly in DECISIONS.md; not surfaced in the tests or ARCHITECTURE 10 |

**10.4 Incorrect or insufficient reference calculations.** Each derived here and checked against the test's own
measurements; none is a defect in NiTCAD's solver, but several are defects in what the tests claim.

- **I-1 `equilibrium:156` treats an exact solution as an approximation, and its measured gap is the metric's bias
  (insufficient).** For a *symmetric* abrupt junction (N_A = N_D = N, Boltzmann, constant ε and n_ie), the first
  integral of Poisson's equation on the n side, with u = (ψ − ψ_bulk)/V_T, is
  E² = (2qV_T/ε)[N(e^u − 1 − u) + p_n(e^{−u} − 1 + u)], p_n = n_ie²/N. By symmetry the junction sits at
  u_j = −V_bi/(2V_T), and since e^{V_bi/2V_T} = N/n_ie, p_n e^{−u_j} = n_ie. Hence
  E_j² = (2qN/ε)(V_bi/2 − V_T + V_T e^{−V_bi/2V_T}) + O(qn_ieV_T/ε)
  = (qN/ε)(V_bi − 2V_T)·(1 + O(n_ie/N) + O(e^{−V_bi/2V_T})), which is the test's "depletion approximation with the
  2 V_T correction" exactly, up to relative terms of 1e-6 (1e16) to 1e-8 (1e18). So the reference is class E for the
  symmetric cases, and its 0.5% tolerance is 3–50 times the measured differences. Those differences (0.01%, 0.10%,
  0.17% at 1e16, 1e17, 1e18) are numerical: the metric, the largest edge difference quotient, is the mean field over
  the junction's first edge, low by (qN/ε)(h/2) with h = 1e-8 cm, i.e. 0.02%, 0.07%, 0.21% of E_j: the same size as
  the measured gaps. The asymmetric case (1e18/1e16, asserted only as "> 1.2× the approximation") is a genuine A-class
  failure of the approximation (the heavy side's majority tail spills a charge comparable to the depletion charge into
  the light side; first-integral estimate here: the field rises well above the depletion value, consistent with the
  measured +39%) and is exactly computable by the same first integral (already coded in `heterojunction:116`).
- **I-2 `bias:165` ideal-diode law: inconsistent assumptions inside a wide band (insufficient).** The reference is the
  short-base Shockley current, but (a) SRH is on (legacy options), so the depletion-region recombination current
  (Sah–Noyce–Shockley), which the formula omits, is in the simulated current; (b) the neutral widths use the
  *zero-bias* depletion width (0.147 µm here) while the comparison is at 0.5 V (0.093 µm), a 2.9% change of J_s; (c) D
  comes from NiTCAD's `caughey_thomas_mobility` (C-9). The ±15% band absorbs all three without measuring any. The
  ideality check (±0.02 from 0.3 V) is sound.
- **I-3 `bias:284` reverse leakage exponent in (0.5, 1.5) (insufficient).** For SNS generation in the depletion region
  J ∝ W ∝ (V_bi − V)^{1/2}, so the expected exponent is 0.5 (plus small terms from the neutral regions' diffusion
  current); the window also admits an exponent of 1, which no mechanism here predicts. V_bi is hard-coded (0.83 V).
- **I-4 `bias:257` two meshes, 3% (insufficient as solution verification).** No order, no Fs = 3 GCI (two-grid
  factor, NPARC), no third mesh.
- **I-5 `small_signal:856` depletion capacitance: the 0.88% is attributed, not shown.** For the 1e19/1e16 junction at
  −2 V the candidate terms are of the right size ((L_D/W)² ≈ 0.45%, the p⁺ side's penetration 0.1%, the heavy side's
  majority tail), but the exact dQ/dV from the first integral (as in I-1) would settle it; tolerance 1.5%.
- **I-6 `mosfet:222` legacy bands.** V_th ∈ (−0.5, 1.5) V and SS ∈ (55, 120) mV/dec are qualitative; the added
  swing bracket (depletion capacitance at φ_s = φ_F and 2φ_F) is a sound A-class interval. No Pao–Sah or charge-sheet
  reference exists.
- **I-7 Order estimates against a finite reference (`unstructured:119`, `unstructured3d:302`, FC8, UC7).** Orders are
  measured against a tensor solve 4× finer than the finest mesh. If that reference converges at first order in H1, its
  own error is about a quarter of the finest unstructured mesh's, which biases the observed order; the bias is not
  estimated, and no asymptotic-range check (NPARC's GCI₂₃ / (r^p GCI₁₂) ≈ 1) is made. An exact (manufactured) solution
  removes the question.
- **I-8 `energy_balance:375` (EC10) compares different models.** Read in Aste and Vahldieck's arXiv text (v1 HTML):
  their energy-balance result of 1.87 v_s uses the Baccarani–Wordeman *temperature-dependent* τ_w (their eqs 15–16) and
  a mobility with μ_min 80, Δμ = 1430 − μ_min, N_ref 1.12e17, exponent 0.72 (their eq 36, consistent with the test's
  80 + 1350/(1 + (N/1.12e17)^0.72)); NiTCAD uses a constant τ_w (documented in ARCHITECTURE 11.3). **The test comment's
  "beta = 2" is not in the source text read** (no β appears; their field dependence enters through τ_p(T)). β = 2 is a
  model input, not only a comment: the fixture sets the Caughey–Thomas field exponent `.beta = 2.0`
  (`energy_balance:125`), and both comments that list it as the source's (`energy_balance:117`, `:377`) attribute
  NiTCAD's parameter choice to the paper. The gate (peak > 1.3 v_s, located in the last 30 nm of the n region) is
  qualitative, which is appropriate for a model mismatch, but neither comment should cite the source for β.
- **I-9 `small_signal:994` Shockley admittance: correct and now quantified.** Shockley (BSTJ 28, 1949, eq 4.24 and §4.1–
  4.2, read in the archive.org OCR) gives (1 + iωτ_p)^{1/2} for the *hole-diffusion term* A_p only, for a long, uniform
  n region, a narrow transition region, low injection and body recombination; the full junction admittance (eq 4.23)
  adds the electron term and the transition capacitance. The test normalizes the real part, which removes C_T, and
  attributes its 3.2e-3 to the electrons injected into the p⁺ side. Estimated here: the electron (short-base, transit
  W_p²/2D_n ≈ 2 ns ≪ τ) share of G₀ is about 7.6e-3 (D_n ≈ 2.6 cm²/s at 1e19, W_p = 1 µm, L_p ≈ 34 µm), which shifts
  the normalized ratio by ≈ 7.6e-3 × (1 − 1/Re√(1+3i)) ≈ 2.4e-3 at ωτ = 3: the attribution holds in size.
- **I-10 `small_signal:735` conductance method: verified.** Nicollian and Goetzberger (BSTJ 46, 1967, §4.3, OCR):
  "G_p/ω goes through a maximum when ωτ = 1 … The value of G_p/ω at the maximum is C_s/2"; their eq (17) itself is
  illegible in the OCR, so the full curve's form is consistent with, not quoted from, the source.

**10.5 Missing theoretical coverage (against §6's catalogue).** Already present in some form, to be catalogued and
re-homed (with the original gates kept, TVD-14): TVB-02/03 (`references.hpp`), TVB-06, 08, 09, 11, 12, 13, 22, 24,
27, 28 (K2 above), 16, 17, 18, 21, 23 (K3, with I-1–I-6), 31 (K4), 33, 34 (K7). **Missing entirely:** manufactured
solutions for Poisson and drift-diffusion on any mesh (TVB-04, 14); curved-geometry exact solutions on unstructured
meshes (TVB-05, 29); the exact pn-junction first integral as the reference for depletion tests (TVB-07: the code exists
for heterojunctions only); Pao–Sah / charge sheet (TVB-20); Haynes–Shockley (TVB-19); exact Airy-function tunnelling
transmission against WKB (TVB-25); Schrödinger–Poisson against MLDA (TVB-30); passivity, Kramers–Kronig and systematic
limit continuity (C8–C10); any GCI with safety factor and asymptotic-range check; any error-interaction measurement
(owner condition 2).

**10.6 Literature verification log (TVR-0).** Read 2026-10-31 by an automated fetch of open copies; quotes kept in
the session record.

| Source | Accessed | Status | Effect |
|---|---|---|---|
| NASA NPARC, "Examining Spatial (Grid) Convergence" | Read (HTML; equations as images, worked forms as text) | **Verified:** Fs = 1.25 for three or more grids and 3.0 for two; the asymptotic check GCI₂₃ ≈ r^p GCI₁₂ (worked ratio 1.002); p from three grids at constant r (worked form) | §4.2's GCI and asymptotic check stand; I-4, I-7 |
| Shockley, BSTJ 28, 435 (1949) | Read (OCR) | **Verified with scope:** (1 + iωτ_p)^{1/2} is the hole-diffusion admittance A_p; assumptions §4.1–4.2 as in I-9; total admittance eq 4.23 | I-9; TVB-18 must use eq 4.23's full form |
| Nicollian and Goetzberger, BSTJ 46, 1055 (1967) | Read (OCR) | **Peak verified** (§4.3); single-level formula eq (17) illegible; continuum peak location not found in the text | I-10; TVB-22 keeps the single-level form, continuum deferred |
| Aste and Vahldieck, arXiv physics/0312021 | Read (HTML v1; PDF not extractable) | **Verified:** device, μ₀(N), v_s, T-dependent τ_w (eqs 15–16), 2.61 v_s (HDM) and 1.87 v_s (EBM) at 1 V; **β not stated**; peak position only in a figure | I-8; TVB-33 |
| Salari and Knupp, SAND2000-1444 | **Not accessible** (OSTI refused connections) | Not verified | §5.2's MMS method stands on the general literature only; to be read before TVR-2 |
| Eymard, Gallouët, Herbin, *Finite Volume Methods* | **Not accessible** (bot challenge) | Not verified | §5.3's orders marked [unverified]; to be read before TVR-2 |
| Masetti 1983; Altermatt 2006; Lombardi 1988; Paasch–Übensee (MLDA); Sze–Gibbons 1966; Pao–Sah 1966; Brews 1978; Gundlach 1966; Simmons 1963; Crowell–Sze 1966 | Not attempted in TVR-0 (paywalled per DECISIONS.md, or needed only in later phases) | — | C-10; each phase reads its sources before its proposal |

**10.7 Prioritised corrections and recommendations.** None changes an existing test, tolerance or ARCHITECTURE.md
without a ruling; tightening is done by *adding* criteria, never by editing a passing gate silently (TVM2).

| Pri | # | Recommendation | Kind | Addresses |
|---|---|---|---|---|
| P0 | R-1 | Relabel, in ARCHITECTURE 10 and the tests' comments, every K5/K6/C-1–C-8 check by what it is (consistency, self-convergence, legacy regression); keep them as pins. ARCHITECTURE 10's "Built-in potential" row already states it does not test the solve; it gains only the n_ie circularity (C-1). The "Diode J(0.5 V)" row already records the analytic 1.321e-2; it gains the label "legacy cross-version regression, not theory verification" (C-7). | Documentation (needs approval: architecture text) | C-1–C-8 |
| P0 | R-2 | Correct this proposal's unread claims: §5.3's Eymard–Gallouët–Herbin orders (done here, marked [unverified]); obtain and read EGH and Salari–Knupp before TVR-2. | Documentation; source acquisition | 10.6 |
| P0 | R-3 | Provenance register: for every M/A reference, the PA–PD grade, the sources actually read, and the transcription path (DECISIONS → test); surface PB/PC grades in test output and ARCHITECTURE 10.1. Ask the owner for the paywalled primaries (Masetti, Altermatt, Lombardi, Paasch–Übensee, Sze–Gibbons, Pao–Sah) where obtainable. | Documentation; source acquisition | C-10, owner condition 3 |
| P1 | R-4 | Exact first-integral references for the pn junction: equilibrium field and potential (symmetric and asymmetric), and dQ/dV at reverse bias; a bias-free field metric (D at the junction from the dual-face flux, or the edge field extrapolated with the known slope ρ/ε); tolerances from U_num on a ladder. Added beside `equilibrium:156` and `small_signal:856`, which stay. | New references (TVR-2) | I-1, I-5 |
| P1 | R-5 | Diode family with consistent assumptions: SRH off for the Shockley comparison (or SNS added), W(V) not W(0), parameters from the benchmark spec; the reverse-leakage exponent against the SNS integral; keep `bias:165`, `:284` as they are. | New references (TVR-4) | I-2, I-3 |
| P1 | R-6 | Manufactured solutions (Poisson, then coupled drift-diffusion) on tensor, 2D, 3D and swept meshes, with the TVQ3 hook; orders with NPARC's GCI and asymptotic check; the K5 gates kept as solution verification. | New references (TVR-2) | C-5, I-7, 10.5 |
| P1 | R-7 | Input independence: reference-side parameter sets written from the sources (or DECISIONS.md's transcriptions, with their grade), never the implementation's structs or derived-quantity functions, for every K2/K3 reference; where sharing is deliberate (code verification of one term), say so in the test. | Reference rework (TVR-1/2) | C-9 |
| P2 | R-8 | EC10: remove the source attribution of β = 2 from both comments (`energy_balance:117`, `:377`), record β = 2 (`:125`) as NiTCAD's parameter choice, state the model mismatch (constant vs T-dependent τ_w), keep the qualitative gate; a reference-only Baccarani–Wordeman τ_w would make a P-class comparison meaningful. | Documentation; optional reference | I-8 |
| P2 | R-9 | Shockley admittance with the electron term and C_T (eq 4.23), so the 3.2e-3 becomes a predicted, not attributed, residual. | Reference (TVR-4) | I-9 |
| P2 | R-10 | Pao–Sah and charge-sheet references for the long-channel MOSFET (sources to be read first). | New reference (TVR-4) | I-6 |
| P2 | R-11 | Consistency suite additions: passivity, Kramers–Kronig (diagnostic, band-limited), limit continuity sweeps; every check with NOT APPLICABLE predicates. | New checks (TVR-1) | 10.5, owner condition 1 |
| P2 | R-12 | Error-interaction ladders (h × tol, and Δt for transients) on the P1 benchmarks. | Method (TVR-2) | Owner condition 2 |

**10.8 Questions for the owner (before TVR-1).**
- TVR0-Q1: Approve R-1's relabelling (a documentation change to ARCHITECTURE 10 and test comments only)?
- TVR0-Q2: Can the paywalled primaries in R-3 be supplied (PDFs placed locally), so the PB/PC references can be raised
  to PA? Without them, those M-class references stay "provenance PB/PC" in every report.
- TVR0-Q3: R-7 changes how existing references obtain their inputs. Apply it to new references only (recommended),
  or also rework existing tests' expected values (each such change would be its own reviewed commit)?
- TVR0-Q4: Proceed to TVR-1 (reference numerics, metrics, report, consistency suite) as the next proposal?

**Owner rulings (2026-11-04):** "TVR0-Q1 to Q4 in §10.8, choose recommended." Read as:
- **Q1 approved:** R-1's relabelling (ARCHITECTURE 10 rows and test comments, documentation only; no test, gate or
  tolerance changes), as amended by the review: the "Built-in potential" row gains only the n_ie circularity, the
  "Diode J(0.5 V)" row the legacy-regression label, and the EC10 comments lose the β attribution (R-8).
- **Q2:** no recommendation was stated; the default holds. The PB/PC references keep their grade in every report
  until the owner supplies the primaries; supplying them stays open.
- **Q3 ruled as recommended:** R-7 (input independence) applies to new references only; existing tests' expected
  values are not reworked.
- **Q4 approved:** TVR-1 (reference numerics, metrics, report, consistency suite) is the next proposal.

**10.9 Applied under Q1 (2026-11-04, uncommitted).** ARCHITECTURE 10's "Built-in potential" row gains the n_ie
circularity (C-1) and the "Diode J(0.5 V)" row the legacy-regression label (C-7). A one- to three-line comment
"TVR-0 label (DECISIONS.md Unit 55 §10, …)" is added at the top of 19 test cases: `equilibrium:92` (C-1);
`bias:144` (C-7), `:257` (I-4); `impact_ionization:170`, `:203` (C-2); `deck/schottky:286` (C-3); `small_signal:232`,
`:271`, `:681`, `:705`, `:783` (C-4); `unstructured:119`, `:200`, `:356`, `unstructured3d:302`,
`complex_convergence:396` (C-5); `transient:209` (C-6); `mosfet:222` (I-6); `mos:172` (C-8). The EC10 comments
(`energy_balance:117`, `:377`) no longer attribute β = 2 to the source (R-8). Comment-only: no test, gate,
tolerance or code line changed (checked on the diff: every added or removed test line is a comment). Not labelled:
`electrothermal:602` (G-PARABOLA compares with an exact parabola; the audit's K6 listing reflects its legacy origin,
not a circularity) and the K5 gates named only in passing (`mosfet:499`, `complex_convergence:452`,
`electrothermal:1334`), left for the catalogue links of TVR-1.

#### 11. TVR-1 (Unit 55a) proposal: reference numerics, metrics, contract, report, consistency suite

**Status:** design for review (2026-11-04). No code until the owner rules TVR1-Q1–Q6 below and says go. Branch
`verify/theoretical-validation` (the proposal), implementation on a new branch per unit.

**11.1 Scope.** In: (a) the `verify/reference` numerics of §3.2 with their L0 tests; (b) the reference-model contract
(§3.3) as types; (c) the metrics of §4.2 as one library; (d) outcomes (§4.4, TVD-16) and the verification report
(TVD-11); (e) `nitcad verify`; (f) `verify/consistency` C1–C12 (§5.1, R-11) with applicability predicates; (g) the
catalogue entries TVB-01, 02, 03 and 32; (h) the catalogue links for the existing tests labelled in TVR-0. Out: mesh
and tolerance ladders and GCI driving (`verify/study`, TVR-2, though the GCI *formulas* are in (c)); the TVQ3 MMS hook
(TVR-2); device benchmarks other than TVB-01–03; the Sturm/inverse-iteration eigensolver (TVQ4 ruling; first needed by
TVB-30 in TVR-5, so built there with its own L0 tests); any change to existing tests, gates, tolerances or the solver.

**11.2 Libraries and dependency rules (TVD-19).**

| Target | Depends on | Never includes | Enforced by |
|---|---|---|---|
| `nitcad_verify_reference` | the standard library; `base`'s constants (except TVB-01, whose reference holds its own CODATA literals) | `physics`, `assemble`, `solve`, `mesh`, `device`, `linalg`, `results`, Eigen | link list + a header-boundary test (`#if __has_include` / `#ifdef` guards, as `tests/linalg/header_boundary_test.cpp`) |
| `nitcad_verify_consistency` | `results`, `mesh` (dual volumes, edges, couplings), `device` (doping, permittivity, contacts) read-only; `verify_reference` (metrics, small eigenproblems) | `physics`, `assemble`, `solve`, `linalg` | the same |
| `nitcad_verify_report` | `verify_reference`, `results` (Provenance), `io` (checksummed records) | `physics`, `assemble`, `solve` | link list |
| `nitcad_verify_benchmark` | everything (it runs NiTCAD: the subject of the check) | — | — |
| `nitcad` (cli) | adds `verify` | — | — |

The consistency layer judges NiTCAD's results using only the results, the mesh and the device: it never re-evaluates
NiTCAD's physics, so a consistency check cannot be satisfied by the same function that produced the result. Where a
check would need physics (n_ie for np = n_ie², the rates' own scales) it is restated in quantities the results carry
(below), or it is NOT EVALUATED.

**11.3 Reference numerics and their L0 criteria (frozen before code).** Each primitive cites its source and is graded
PA–PD (§3.3). Tolerances are in units of u_dd = 2⁻¹⁰⁴ ≈ 4.9e-32 (double-double unit roundoff) or ε = 2⁻⁵² for double.

| # | Primitive | Source (to be read in full before code) | L0 criteria |
|---|---|---|---|
| TC1 | Double-double arithmetic: two_sum, two_prod (FMA), +, −, ×, ÷, sqrt; exp, expm1, log, pow | Dekker (Numer. Math. 18, 1971); Hida, Li and Bailey (QD library report, 2000) | Constants π, e, ln 2, √2 to 32 digits against literals from OEIS A000796, A001113, A002162, A002193 (PA); identities exp(log x) = x, exp(a + b) = exp a · exp b, (√x)² = x, (x / y) · y = x to 8 u_dd relative on 10⁴ log-spaced points over each function's range; error-free transformations exact (hi + lo reproduces the exact sum/product, checked by an integer-exact construction) |
| TC2 | Cross-check with the existing test helper | `tests/physics/reference_arithmetic.hpp` (TVR-0: correct as written) | The new library and the existing header agree to 2 u_dd on the TC1 grid. The header is not changed or included by the new library (TVD-14); the agreement test lives in the verify tests. (TVR1-Q1) |
| TC3 | Gauss–Kronrod 7/15, adaptive, with error estimate | Piessens et al., QUADPACK (1983); the qk15 nodes and weights from netlib's public QUADPACK source (PA) | Nodes and weights to their printed digits; monomials xᵏ on [−1, 1] exact to 4 u_dd for k ≤ 22 (the 15-point Kronrod rule's degree 3n + 1) and not for k = 24 (shows the degree is real); adaptive: on ten integrals with closed forms (smooth, peaked, oscillatory), the true error ≤ the reported estimate in every case, and ≤ the requested tolerance |
| TC4 | tanh–sinh (double-exponential) quadrature | Takahasi and Mori (1974); Bailey, Jeyabalan and Li (2005) | ∫₀¹ x^{−1/2} = 2, ∫₀¹ ln x = −1, ∫₀¹ x^{−1/2} ln x = −4, ∫₀^∞ e^{−x²} = √π/2 to 16 u_dd; the error estimate bounds the true error on each |
| TC5 | Bracketed Brent root finding | Brent, *Algorithms for Minimization without Derivatives* (1973), ch. 4 | Roots of x³ − 2x − 5 (Newton's cubic, root known to 40 digits, PA from the literature), cos x − x, and e^x − 1 near 0 to 4 u_dd; refuses an unbracketed interval; iteration count within Brent's bound |
| TC6 | Dormand–Prince 5(4) with dense output | Dormand and Prince (J. Comput. Appl. Math. 6, 1980); Hairer, Nørsett and Wanner, vol. I | Tableau to printed digits; observed order 5.0 ± 0.2 on y′ = y and the harmonic oscillator at fixed steps; adaptive global error proportional to the tolerance over four decades; dense output order 4 |
| TC7 | Chebyshev collocation with domain splitting | Trefethen, *Spectral Methods in MATLAB* (2000); Boyd, *Chebyshev and Fourier Spectral Methods* (2001) | −(εu′)′ = f on two subdomains with an ε jump and an exact piecewise solution: error decays geometrically (each doubling of N reduces it by ≥ 10³ until it reaches 64 u_dd·max\|u\|); interface conditions (continuity of u and εu′) hold to the same level |
| TC8 | Complex-step differentiation (double) | Squire and Trapp (SIAM Rev. 40, 1998); Martins, Sturdza and Alonso (2003) | Derivatives of e^x, sin x, x^{3/2}, and a nested composite to 4ε relative at h = 1e-200; agreement with the analytic derivative where one exists |
| TC9 | Small Hermitian eigenproblems (n ≤ 16): cyclic Jacobi | Golub and Van Loan, *Matrix Computations* (4th ed.), §8.5 | Known spectra (diagonal, rank-one updates, the Hilbert matrix's smallest eigenvalue to its published digits); residual ‖Av − λv‖ ≤ 8nε‖A‖; used only by C8 (passivity) |

Each primitive reports an error estimate, which the contract carries as U_ref. **Mutation requirement (house
practice):** each L0 test set must fail on at least one deliberate corruption of its primitive (a wrong weight, a
dropped lo term, a wrong tableau entry), recorded in the as-built.

**11.4 Contract and outcome types.** `verify::ReferenceInfo` (id, class E/D/A/M/H/P/X, independence grade I0–I3,
provenance grade PA–PD, citation with equation numbers, version, the G1–G15 subset, assumptions, applicability
predicates); `verify::Predicate` (name, evaluated value, threshold, holds); `verify::Quantity` (name, unit, value,
U_ref); `verify::Outcome` = pass | finding | fail | not_applicable | not_evaluated, with the predicate or missing input
named for the last two (TVD-16). A reference graded PC cannot produce `pass` (TVD-18): its comparisons report
`finding` with "provenance insufficient". These are plain data types with no behaviour beyond validation of their own
fields.

**11.5 Metrics library (TVD-7).** The formulas of §4.2: relative error with floor; log error with floors; discrete
L∞, dual-volume L2 and dual-face H1 norms (from `mesh` geometry, in the consistency/benchmark layers; the reference
layer only gets vectors and weights); observed order from three errors, or from three solutions at constant r;
asymptotic-range check; GCI with F_s = 1.25 (≥ 3 grids) and 3.0 (2 grids); the factorial interaction contrast ε_int
of §4.1 and the combination U = Σ\|U_k\| + \|ε_int\|. **L0 criteria:** exact recovery of p from synthetic e_k = C hᵏᵖ
(to 1e-12); the NPARC page's worked example reproduced to its printed digits (verified source, §10.6); ε_int = 0 for
an additive synthetic e(h, tol) and the exact cross-term for a multiplicative one; GCI refused (`not_evaluated`) when
the three solutions are not monotone or p̂ ≤ 0.

**11.6 Consistency suite (`verify/consistency`).** Each check is a function of (result, mesh, device) returning an
Outcome, the measured residual, the tolerance and the predicates. Tolerances use quantities the result already
carries; none is set after a run.

| # | Restated in result quantities | Predicates (else NOT APPLICABLE) | Tolerance (frozen) |
|---|---|---|---|
| C1 | Terminal and edge currents zero; quasi-Fermi levels flat and equal (F_n = F_p everywhere, which is np = n_ie² under the run's own statistics, with no physics call) | all biases 0; no optical or prescribed generation; no transient | currents ≤ their `terminal_current_resolution`; edges ≤ that resolution; \|F_n − F_p\| ≤ 64ε·max(\|E_c\|, \|E_v\|, V_T) per node |
| C2 | Each stored mechanism's rate integrates to nothing measurable: q∫\|R_k\| dV ≤ the smallest terminal resolution | C1's predicates; rates stored (else NOT EVALUATED) | as stated (a rate below it cannot drive a resolved current) |
| C3 | Dual-cell Gauss law: Σ_faces ε_f (A_f/d_f)(ψ_i − ψ_j) = the cell's charge q(p − n + N_D − N_A)V_i plus interface and trap charge; global neutrality with electrode charges | converged state; Boltzmann or FD doping fully ionized, or ionized densities available (else NOT EVALUATED for incomplete ionization, which the results do not store) | the Newton tolerance scaled as the solver states it, plus 64ε times the largest face flux |
| C4 | Σ terminal currents = 0 (DC); with displacement current (transient; AC column sums) | — | Σ of the contacts' resolutions (DC, transient); `resolution[k][j]` (AC) |
| C5 | Σ I V = boundary heat + stored energy rate (electrothermal) | electrothermal on; steady; `thermal_contact_heat` present | the thermal residual's stated scale; NOT EVALUATED where the volume heat terms are not in the results (listed gap) |
| C6 | Gate charge: ∫ I_g dt = ΔQ_g by the run's own quadrature order | transient with a gate | the step's local error estimate summed |
| C7 | \|Y_ij − Y_ji\| / max\|Y\| | all biases 0; no magnetic field (none in NiTCAD); AC | (resolution_i + resolution_j) / max\|Y\| |
| C8 | λ_min((Y + Yᴴ)/2) ≥ −tol, by TC9 | all biases 0 | ‖resolution column‖ |
| C9 | Re Y(ω_min) = dI/dV and Im Y(ω_min)/ω_min = dQ/dV from a pair of DC points the caller supplies; band-limited Kramers–Kronig on a log-spaced Y(ω) | linear regime (the caller states the amplitude); for KK at least 8 points per decade over ≥ 4 decades | DC part: the DC points' resolution plus the finite-difference error (stated); **KK: diagnostic only**, reported as a measured value with the tail-truncation estimate, never pass/fail |
| C10 | Difference to the limit model decays at the stated rate over a caller-supplied sequence | the limit parameter is swept (≥ 3 values) | rate within 0.2 of the stated one |
| C11 | Mirrored or rigidly moved devices: the caller supplies the node map; terminal quantities equal | the map is a symmetry of the device (checked on doping and geometry) | 64ε relative |
| C12 | n, p > 0; T > 0; occupancies in [0, 1] | always | exact |

**Each check must be shown to detect what it claims:** a defect-injection test per check (a corrupted current, a
shifted F_n on one node, a broken Y symmetry, a negative density, a dropped interface charge) fails it, and real solves
of three existing decks (`examples/pn_diode_1d.deck` (DC sweep), `examples/moscap_cv.deck` (small-signal C–V, two
contacts), `decks/schottky_diode_1d.deck` (DC sweep)) pass it or report NOT APPLICABLE with the right predicate. C7
and C8 need an AC point at zero bias: where the C–V sweep has none, the benchmark adds one through the API (the deck
is not edited). These real-solve runs are the TVB-32 entries of TVR-1; later phases add their
benchmark devices.

**11.7 Report and CLI (TVD-11).** The report is written in `io`'s checksummed record conventions (Unit 42
provenance: build, commit, dirty flag, platform, thread count) and a CSV summary: per benchmark, the reference's
identity and grades, each quantity with q_sim, q_ref, U_ref, the U_num components and ε_int (or "not estimated"), the
metric, the criterion, the outcome and its reason. Two runs on one machine and thread count give bit-identical reports
except the wall time (the same rule as Unit 42's digest). `nitcad verify [--tier fast|full|nightly] [--id TVB-…]
[--report <file>]` follows the runner's existing convention (`src/cli/main.cpp`): 0 when no benchmark fails, 1 on any
`fail` or a run or I/O failure (the reason on stderr), 64 on usage; `finding`,
`not_applicable` and `not_evaluated` do not change the exit code but are counted in the summary line, and an
all-`not_applicable` benchmark is printed as "unverified". (TVR1-Q4.)

**11.8 Catalogue entries in TVR-1.** TVB-01 (constants, bit-equal to CODATA 2018 decimal strings typed from NIST's
table, PA); TVB-02 (F_{1/2}, F_{−1/2}, F_{3/2} and the inverse against TC3/TC4 double-double quadrature, ≤ 1e-12
relative over η ∈ [−40, 40], with Blakemore's approximation as a printed-digits cross-check); TVB-03 (Bernoulli
B(x) ≤ 4ε relative); TVB-32 as above. The existing double-double gates (`physics/fermi_dirac:69`, `assemble/flux:27`)
stay as they are; TVB-02/03 are new and their agreement with the existing ones is recorded, not asserted (TVD-14).
The TVR-0 labels are linked from the catalogue as `existing: <file:line> (kind, class, grade)`.

**11.9 Tiers and cost (TVD-13).** All of TVR-1 is `fast` (Debug-safe): the L0 tests are scalar, the consistency
checks are linear passes over results, and the three deck solves are 1D. Budget ≤ 60 s added in Debug; measured at
build and recorded. A new ctest `verify` runs in both configurations.

**11.10 Acceptance criteria (VR1–VR10, frozen for review).**
- VR1: TC1–TC9 pass their L0 criteria, each with its recorded mutation.
- VR2: the header-boundary tests for `verify_reference` and `verify_consistency` fail to compile when a forbidden
  header is included (shown once by a deliberate include, recorded) and pass otherwise.
- VR3: every reference in TVR-1 carries a filled contract; a PC-graded test reference cannot pass (a test proves it).
- VR4: the metrics' L0 criteria (§11.5) pass, including the NPARC worked example.
- VR5: each of C1–C12 fails on its injected defect and passes or is NOT APPLICABLE (with the right predicate) on the
  three decks; C2, C3, C5 report NOT EVALUATED where the input is missing, not pass.
- VR6: TVB-01–03 pass at their frozen tolerances; any miss is a finding for the owner, and no tolerance moves (TVM2).
- VR7: two report runs bit-identical except wall time.
- VR8: `nitcad verify` exit codes as §11.7, tested.
- VR9: no existing test, gate, tolerance or source line outside the new targets, `CMakeLists.txt` and `src/cli`
  changes (checked on the diff).
- VR10: added Debug time ≤ 60 s; Release and Debug suites green.

**11.11 Design choices.**

| # | Choice | Why |
|---|---|---|
| TVD-19 | The dependency table of §11.2, with header-boundary tests for the two independent layers. | Independence enforced, not promised |
| TVD-20 | Consistency checks are restated in quantities the results carry; where impossible they are NOT EVALUATED, never computed by calling NiTCAD's physics. | A check that calls the checked code is I0 |
| TVD-21 | Kramers–Kronig is diagnostic only. | Band-limited data make a KK pass/fail a tuning exercise |
| TVD-22 | The eigensolver of TVQ4 is built in TVR-5 with TVB-30; TVR-1 has only the n ≤ 16 Jacobi for C8. | Build each primitive where it is first used and tested |
| TVD-23 | `finding`, `not_applicable` and `not_evaluated` do not fail the CLI; `fail` does. | A finding needs an owner ruling, not a red build |

**11.12 Questions for the owner.**

| # | Question | Recommendation |
|---|---|---|
| TVR1-Q1 | Double-double: a new library in `verify/reference`, with the existing test helper left untouched and cross-checked (TC2), or move the existing tests onto the new library? | New library + cross-check; deduplicate later only by ruling |
| TVR1-Q2 | Source acquisition before code: QUADPACK (netlib), Dormand–Prince 1980, Brent 1973, Takahasi–Mori 1974, Hida–Li–Bailey, Golub–Van Loan §8.5. Open copies will be fetched and read; may the paywalled ones (Dormand–Prince, Takahasi–Mori) be supplied, or graded PB from two reproductions? | Fetch open copies; PB where two agreeing reproductions are read, the grade printed |
| TVR1-Q3 | C2's tolerance as the integrated bound q∫\|R_k\| dV ≤ the smallest terminal resolution (computable now), rather than a per-node bound needing a new per-mechanism scale in `results` (a change to existing code)? | Integrated bound now; per-node scale only if a later unit needs it |
| TVR1-Q4 | CLI exit codes as §11.7 (fail → 1, as the runner's failures; finding, NOT APPLICABLE, NOT EVALUATED → 0, counted; usage 64)? | Yes |
| TVR1-Q5 | The three decks for TVB-32 (`examples/pn_diode_1d`, `examples/moscap_cv`, `decks/schottky_diode_1d`), with a zero-bias AC point added through the API for C7/C8? | Yes; TVR-2 adds its benchmark devices |
| TVR1-Q6 | Branch `verify/reference-numerics` from main after the TVR-0 and §11 documentation is committed, one commit per TC group, local only until you authorise a push? | Yes |

**Owner rulings (2026-11-04):** "choose recommended". Q1: a new double-double library in `verify/reference`, the
existing helper untouched and cross-checked (TC2); deduplication only by later ruling. Q2: open copies fetched and
read; a paywalled source is graded PB where two agreeing reproductions are read, the grade printed. Q3: C2's
integrated bound; no change to `results`. Q4: exit codes as §11.7. Q5: the three decks, with a zero-bias AC point
added through the API. Q6: the TVR-0 and §11 documentation committed locally on `verify/theoretical-validation`;
implementation on `verify/reference-numerics`, one commit per TC group, local only until a push is authorised. Since
the documentation is not on main until this branch is merged, the implementation branch starts from this branch's
commit (so the design it implements travels with it). Implementation starts when the owner says go.

**11.13 Amendments to the frozen criteria (owner rulings 2026-11-04, before any measurement).** The owner
said go on 2026-11-04 ("start TVR-1"). Two TC1 criteria could not be met by any correct double-double
implementation, so they went to the owner before a test ran:
- **TC1 exp identities.** exp(log x) = x and exp(a + b) = exp a · exp b are ill-conditioned: a double-double
  y = log x carries only about 2⁻¹⁰⁶\|y\| absolute precision, and exp multiplies it by \|y\| (≈ 170 u_dd at
  x = 1e300 for a correctly rounded log). Ruling: the tolerance is 8 u_dd × the identity's condition number,
  max(1, \|ln x\|) and max(1, \|a\| + \|b\|), over the full range. (√x)² = x and (x / y) · y = x stay at 8 u_dd.
- **TC1 constants "to 32 digits".** A double-double holds about 31.9 digits; half a unit in the 32nd digit of
  ln 2 is 7.2e-33 relative, below what a correctly rounded double-double can always reach. Ruling: \|c − L\| ≤
  1 u_dd · \|L\| against the 40-digit OEIS literal, compared exactly with big integers.
- **Domain (recorded, not a ruling).** Double-double precision exists only where the lo part is a normal
  double; the TC1 ranges are \|x\| ∈ [1e-280, 1e280] (quotients and intermediates too). The first run used
  1e-300 and failed there by subnormal lo parts, not by the arithmetic; the header's stated domain was wrong
  and was corrected.
- **TC2 log and pow (after measurement; two rulings).** The first TC2 run missed 2 u_dd on log (2.9 u_dd at
  x = 0.94) and pow (128 u_dd where \|b ln a\| ≈ 267). The log difference was attributed to the helper: on
  expm1(log(1 + d)) = d (\|d\| ∈ [1e-30, 0.1], conditioned by 1) the new log's worst error is 1.41 u_dd and the
  helper's 8e15 u_dd at d = −5.8e-17. Ruling 1: log is compared absolutely, \|Δ\| ≤ 2 u_dd · max(1, \|ln x\|) (the
  helper's guarantee), pow relative to 2 u_dd · max(1, \|b ln a\|). The rerun still missed pow (8.3 per unit at
  a = 1.024): the helper's absolute log error is multiplied by \|b\| inside exp(b log a). Ruling 2: pow relative
  to 2 u_dd · max(1, \|b\|, \|b ln a\|).
- **New TVR-0 finding I-11 (the helper, `tests/physics/reference_arithmetic.hpp`, unchanged per TVD-14).** Its
  log is Newton on exp from std::log(hi) and is accurate to ≈ u_dd *absolutely*; relative to its own value it
  degrades near x = 1, to double precision for \|x − 1\| < 2⁻⁵³ (8e15 u_dd measured). Its pow inherits \|b\| times
  that. Whether any existing gate evaluates the helper's log or pow near 1 at a relative tolerance is not yet
  audited; open for a ruling (no existing test or gate changed).
- **TC3 "nodes and weights to their printed digits" (after measurement).** dqk15.f prints 33 digits, beyond a
  double-double's. Measured: nodes ≤ 0.14 u_dd; Gauss weights ≤ 0.95 u_dd (after replacing 1 − x² by
  (1 − x)(1 + x), which had cost wg(1) 1.76 u_dd near x = 0.949); Kronrod weights 0.04–2.41 u_dd, from an 8 × 8
  moment solve. Ruling: nodes 1 u_dd (as the TC1 constants), weights 4 u_dd (TC3's own monomial unit).
- **TC4 secondary epsilon (after measurement).** With BJL's usual ε₂ = ε₁² (weights to ≈ 1e-63) the endpoint term
  d4 of the published estimate stays at ≈ 7e-29 for ∫₀¹ x^{−1/2} ln x, above its 16 u_dd tolerance (3.2e-30), so
  the run ended `limit_reached` with a true error of 4.7 u_dd. Bailey's 2006 note allows ε₂ = ε₁ⁿ for any n ≥ 2.
  Ruling: n = 3 (weights to ≈ 1e-94, abscissa distances to ≈ 1e-96, inside the double-double domain).
- **TC4 "the error estimate bounds the true error on each" (after measurement).** With n = 3 the published
  estimate fell below the true error twice: ln x 4.36e-32 against 4.37e-32, x^{−1/2} ln x 1.36e-31 against
  2.89e-31. Its floor d3 = ε₁ · max\|w f\| assumes f exact to ε₁; here log and sqrt are themselves accurate to a few
  u_dd. Ruling: keep the published 10^d unchanged and add one documented term, E = max(10^d, 8 u_dd · h Σ\|w f\|)
  (TC1's 8 u_dd guarantee times the computed ∫\|f\|) for n ≥ 3; a deviation from the sources, recorded in the
  header. Consequence: no relative tolerance below ≈ 8 u_dd is reachable.
- **TC5 "e^x − 1 near 0 to 4 u_dd" (reading, not a ruling).** The root is 0, where a relative error is undefined;
  read as 4 u_dd absolute (f′(0) = 1 sets the scale). All three roots use t = u_dd, for which Brent's (2.18),
  6ε\|ζ\| + 2t, guarantees 2 u_dd at 0; the other two roots are judged relative, as frozen. Open to the owner if
  another reading was meant.
- **TC6 readings (owner rulings before measurement).** "Global error proportional to the tolerance over four
  decades": the least-squares slope of log err against log tol (rtol = atol = tol, five tolerances a decade apart)
  in [0.96, 1.04], the frozen order band 5.0 ± 0.2 divided by 5 (h ∝ tol^{1/5}). "Dense output order 4": the
  continuous extension's order as HNW define it, one step from exact data, the error at θ ∈ {¼, ½, ¾} falling as
  h^{5.0 ± 0.2} (a global dense error reads ≈ 5, since the step values are 5th order; it is printed, not judged).
  Addition approved: the nonlinear y′ = y², y(0) = ½ (linear problems see only the tall-tree conditions), under the
  same bands, and the 17 order conditions of order ≤ 5 checked on the tableau.
- **TC6 ladders (owner rulings after measurement).** With the ladders chosen before measuring (fixed N = 16 … 128;
  tol = 1e-10 … 1e-14): the oscillator's coarsest pair read 4.794 (h = 0.39), y′ = y² read 1.46 / 4.24 / 4.70
  (preasymptotic: 4.86 / 4.94 / 4.97 at N = 128 … 1024) and its slope 0.918 (0.99 over 1e-14 … 1e-18); y′ = y passed
  both. Ruling: one ladder for all three problems, N = 128, 256, 512, 1024, and tol = 1e-14 … 1e-18; the band is
  unchanged.
- **TC6 dense output (owner ruling after measurement).** The two reproductions print different interpolants:
  dopri5.f's (degree 4 in θ, attributed there to Shampine 1986) and Boost.Odeint's (degree 5, cited there to HNW
  pp. 191–192). Both satisfy the eight order-4 conditions to 0.28 u_dd; they coincide at θ = 0, ½, 1 and differ
  elsewhere by ≈ 1e-4·θ²(1 − θ)²(θ − ½) per weight. A test written to assert their agreement was wrong and failed
  (1e22 u_dd); it was replaced. Ruling: dopri5.f's interpolant, graded PB on dopri5.f plus the tests' independent
  derivation (the eight conditions for all θ, b(1) = b, b′(0) = e₁, b′(1) = e₇); the Boost form is asserted to be a
  different order-4 extension equal at θ = 0, ½, 1. Recorded as a finding about the sources.
- **TC7 readings (owner rulings before measurement).** The error is the fine-grid interpolant's (201 points per
  subdomain, Clenshaw on the Chebyshev coefficients), the nodal error printed. The decay rule: on N = 4, 8, 16,
  32, 64 per subdomain, each doubling reduces the error by ≥ 10³ or brings it to the floor; once there it stays;
  the floor is reached by N = 64. Interfaces: the discrete mismatches of u and εu′ at the finest N, and the errors
  of u(c) and of each side's εu′(c) against exact under the decay rule. The problem, fixed before measuring:
  [−1, 0] ∪ [0, 1], ε = 1 and 10, u = eˣ and cos 3x + sin(3x)/30 (u(0) = εu′(0) = 1 on both sides).
- **TC7 floors (owner rulings after measurement).** No correct second-order collocation reaches 64 u_dd·max\|u\| ≈
  3e-30: the error fell geometrically (5.9e-2, 1.2e-5, 2.5e-15, 8.5e-27 at N = 4 … 24) onto a rounding plateau of
  ≈ 10³ u_dd at N = 28–32, 6.8e3 at 64 and 3.9e4 at 128, with the coefficient tail at 1e-32 (Boyd §7.7: condition
  O(N⁴); eq. (7.31)). The interface flux, formed from nodal values, loses a further ≈ N² (2e-28 to 1.4e-27 against
  64 u_dd·max\|εu′\| ≈ 9.5e-29). Ruling: the floors are 64 u_dd·max\|u\|·N² for u and u(c), and 64 u_dd·max\|εu′\|·N²
  for the flux mismatch and the flux errors; the u-continuity mismatch stays at 64 u_dd·max\|u\|. Also approved:
  two of the implementer's own added tests had tolerances (8 u_dd) that ignored the N-term summation and the k
  products forming xᵏ; they are (N + 8) u_dd·Σ\|D_ij\| (Higham's γ_N) for the (6.6) diagonal against Theorem 7
  (measured 9.3 u_dd·Σ against 8 before) and (N + 2k + 8) u_dd·Σ\|D_ij x_jᵏ\| for D·xᵏ (23.9 before).
- **TC8 readings (owner rulings before measurement, 2026-11-04).** (1) The reference is the analytic derivative
  evaluated in double-double (TC1's exp, sqrt; sin and cos by a test-local series after reduction by the library's
  π, checked by sin² + cos² = 1 and a double-angle identity, as TC5), so 4ε judges the complex step and not the
  reference. (2) Points: the composite is Lyness and Moler's f = eˣ/√(sin³x + cos³x) (Martins, Kroo and Alonso
  2000, (15)) at its published point x = 1.5, and 10³ points each of eˣ on [−200, 200], sin x on [−2π, 2π], x^{3/2}
  log-spaced on [1e-50, 1e50] and the composite on [0, 1.5] (ranges where f′·h stays a normal double); 4ε relative
  is judged where \|f′\| ≥ 0.01·max\|f′\| over the grid and printed elsewhere (a relative error is unbounded at a
  zero of f′, and the double input's own position error dominates there). (3) U_ref is Richardson's truncation
  estimate \|d(2h) − d(h)\|/3 (from (14), d(h) = f′ − h²f‴/6 + O(h⁴)), labelled "truncation only; the rounding of
  f is not estimated"; it is tested at h = 1e-3, where it must bound the true error, and at h = 1e-200, where it
  must be at rounding level (≤ 4ε\|d\|). (4) Scope: `complex_step(f, x, h)` plus the functions Table 5 marks as
  new definitions (abs, max, min, sign, dim; choices on real parts), refusing h ≤ 0, non-finite or subnormal h and
  non-finite x. Choice recorded before measuring: x^{3/2} is evaluated as z·√z, not pow(z, 1.5), because a pow
  formed as exp(1.5 log z) carries the condition number \|1.5 ln x\| (TC1), up to 170 on this range; pow is printed.
- **TC8 (owner rulings after measurement, 2026-11-04).** (1) The 0.01·max\|f′\| rule judged only 12 of 1 000 points
  of eˣ and 40 of x^{3/2} (f′ spans ≈ 170 decades there; the implementer's own count check caught it), with every
  point within 1.58ε. Ruling: functions whose f′ has no zero on the range (eˣ, x^{3/2}) are judged at every point;
  the 0.01·max rule stays for sin and the composite; 4ε unchanged. (2) At h = 1e-3 the leading-order Richardson
  estimate fell below the true error by O(h²) relative (eˣ 2.0e-7, x^{3/2} 1.8e-7, composite 3.9e-6; sin above), as
  (14)'s next term predicts when f‴ and f⁽⁵⁾ share a sign (for eˣ: estimate h²/6 − 5h⁴/120, true h²/6 − h⁴/120).
  Ruling: the primitive reports 1.25·\|d(2h) − d(h)\|/3 (F_s of §11.5), and the test adds estimate ≤ 1.25·(1 +
  1e-3)·true so the estimate cannot be arbitrarily loose.
- **TC9 readings (owner rulings before measurement, 2026-11-04).** Golub and Van Loan (a book) and Fettis and
  Caslin (Math. Comp. 21, 1967, the printed Hilbert eigenvalues; AMS 403) were not available. (1) The Hilbert
  reference is OEIS A076823's exact integer coefficients of det(H_n)⁻¹·det(xI − H_n) (PA data), evaluated in
  double-double, with λ_min and λ_max found by TC5's Brent: an independent derivation (TVD-18) in place of
  "published digits". (2) Eigenvalues are judged absolutely, \|λ̂ − λ\| ≤ 8nε‖A‖₂ (the frozen residual's unit;
  Weyl's bound for a backward-stable method), since no double method reaches λ_min(H_n) to its digits; H_n is
  integer-scaled (L·H_n, L = lcm(1, …, 2n − 1), so the input is exact) for n = 3 … 10, and λ_min's relative error
  is printed. (3) Scope: complex Hermitian (std::complex<double>), n ≤ 16, Forsythe and Henrici's rotation (6) with
  (43)–(44) and \|φ\| ≤ π/4, cyclic by rows (III_r) until off(A) ≤ ε‖A‖_F, at most 50 sweeps (else not_converged);
  refuses non-Hermitian input and n > 16; ascending eigenvalues with eigenvectors. Added tests: ‖VᴴV − I‖ ≤ 8nε,
  FH's (42), rank-one updates D + ρzzᴴ against the secular equation solved by Brent in double-double, a complex
  unitary similarity of a known diagonal.
- **Metrics readings (owner rulings before code, 2026-11-04).** Source re-read: NPARC "Examining Spatial (Grid)
  Convergence" (grc.nasa.gov, updated 10 Feb 2021) with its VERIFY run `prD.do`/`prD.out`. Recomputed in double,
  the worked example gives p̂ = 1.7861696, f_{h=0} = 0.9713003, GCI₁₂ = 0.1030826 %, GCI₂₃ = 0.3562493 %, all
  agreeing with the page's printed digits; the page's asymptotic ratio 1.002019 does not (exact 1.0020237): it
  is the quotient of the page's own rounded intermediates (0.356249 / (3.44898 × 0.103083) = 1.0020189). VERIFY's
  output (p = 1.78618479, ratio 0.997980) disagrees with the exact values in the fifth digit, consistent with
  single precision. (1) The L0 test holds p̂, f_{h=0} (and its increment 0.00080), GCI₁₂ and GCI₂₃ to half a unit
  of the page's last printed digit; the ratio is checked twice: from the page's rounded intermediates it must give
  1.002019, and the exact value 1.0020237 is recorded as the page's rounding. VERIFY's output is logged as single
  precision and not used. (2) **Finding against §4.2 / NPARC:** with exactly three grids and p̂ from those grids,
  GCI₂₃ / (r^p̂ GCI₁₂) = f₁/f₂ identically (because (f₂ − f₃)/(f₁ − f₂) = r^p̂), so it tests only that the two
  finest solutions are close. Ruling: the asymptotic-range check needs ≥ 4 grids: successive p̂ within 0.1, and,
  as a second test, the GCI ratio within 10 % of 1; with three grids it is `not_evaluated` ("needs 4 grids"), the
  three-grid ratio computed, reported and labelled "identically f₁/f₂". (3) p̂ from solutions only at constant r
  (ratios equal to 1e-12 relative), else `not_evaluated` ("non-constant r"); the non-constant form (Celik et al.,
  ASME JFE 2008) waits until its source is read and graded.
- **Contract readings (owner rulings before code, 2026-11-04).** (1) Every comparison against a PC reference,
  criterion met or missed, reports `finding` ("provenance insufficient") with its measured values: a PC source can
  neither pass nor fail the code. (2) The types also enforce §2.1 and §3.3: an I0 reference cannot pass (`finding`,
  "regression pin, not evidence"); a contract is invalid when an E reference is not PA or PD with two methods, an
  E reference graded I1 names no second I1/I2 cross-check, an A reference states no known approximation error, or
  a required field is empty; a comparison against an invalid contract is `not_evaluated` with the problem named.
- **C1 (owner ruling after measurement, 2026-11-04).** On `examples/pn_diode_1d.deck` at 0 V the drift-diffusion
  point's \|F_n − F_p\| was 1.4695e-14 eV at node 112 against the frozen 64ε·max(\|E_c\|, \|E_v\|, V_T) =
  1.4395e-14 eV (2 % over); the equilibrium solve's own result was at 3.5e-4 of the bound. The point had taken one
  Newton iteration with a last scaled update of 5.67e-13, and V_T·5.67e-13 = 1.47e-14 eV: the split is the Newton
  stopping error, which C1's frozen tolerance omitted and C3's includes. Ruling: at a drift-diffusion point the
  per-node tolerance is 64ε·max(\|E_c\|, \|E_v\|, V_T) + V_T·u_last, u_last the result's last Newton update (a
  quantity it carries; NOT EVALUATED without a Newton record); the equilibrium solve's result is judged on the
  quasi-Fermi part alone at 64ε (u_last = 0). Recorded as a change after measurement.
- **C3 tolerance (owner ruling before code, 2026-11-04).** "64ε times the largest face flux" leaves out the charge
  term q V (p − n + N_D − N_A), whose operands cancel in a neutral cell, where the face fluxes are near 0 (C1's
  lesson). Ruling: a cell's tolerance is 64ε s_i + N_i u_last, with s_i = Σ_f C_f (\|ψ_i\| + \|ψ_j\|) +
  q V_i (n_i + p_i + N_D,i + N_A,i) (Unit 50's row scale Σ_j \|J_ij x_j\|, plus the doping) and the Newton term
  N_i = Σ_f C_f · 2 V_T + q V_i (n_i + p_i), from the result's last scaled update u_last (ψ in V_T, densities
  relative), as in C1's ruling; C_f = ε₀ ε_f A_f / d_f, ε_f the harmonic mean of the end permittivities (the
  discretization's own). NOT EVALUATED without a Newton record. A reading made before measurement.
- **C3 cells not restated alone (owner ruling before code, 2026-11-04).** Ohmic and Schottky contact nodes are
  Dirichlet rows; their charge is the contact's and is not stored: they are skipped. Each lumped gate and each
  meshed electrode is checked as a group: Σ over its cells, plus its stored `gate_charge` (and, for a lumped gate,
  q Q_f × its patch area), is 0. Each charged interface (fixed charge, traps or surface recombination) likewise:
  Σ over the nodes of its edges (whose interface potential is not stored), plus q σ_f × its area plus its stored
  `interface_trap_charge`, is 0. Groups that share a node are one group; a group containing a contact node joins
  that contact. Global neutrality: Σ space charge + every stored electrode, gate, fixed and trapped charge + each
  ohmic or Schottky contact's charge taken as its cells' Gauss flux is 0, which makes it a sum of the other parts'
  residuals; it is computed separately and labelled that way.
- **C5 tolerance (owner ruling before code, 2026-11-04).** The results store no thermal residual of their own
  (`ConvergenceRecord` mixes every row), and "the thermal residual's stated scale" has more than one reading.
  Ruling: Unit 23's U8 scale restated in result quantities, \|Σ_c V_c I_c − Σ_k H_k\| ≤ E_max Σ_c R_c +
  64ε (Σ_c \|V_c I_c\| + Σ_k \|H_k\|), H_k = `thermal_contact_heat`, R_c = `terminal_current_resolution`, E_max =
  max(\|V_c\|; \|E_c\|, \|E_v\| on semiconductor nodes) + max \|r + 5/2\| k T_max, r the thermopower exponents of
  the device's semiconductors (copied into `DeviceData`, a material constant like the permittivity) and T_max the
  largest stored lattice temperature. No separate Newton term: U8 measured the Newton tolerance well inside E_max R_I.
  (U8 itself used 3kT, which assumes r + 5/2 ≤ 3.) A reading made before measurement.
- **C5 volume heat (owner ruling before code, 2026-11-04).** The per-node heat ∫H dV is not in the results;
  if it were a part of the check, C5 would always be NOT EVALUATED. Ruling: C5 judges the global balance only; the
  local part is a listed gap in the contract (an assumption) and in the as-built notes, not a part.
- **C5 scope (owner ruling before code, 2026-11-04).** Steady `BiasPoint`s only (a transient is outside "steady");
  predicates "electrothermal on" (the run's model, in `DeviceData`) and "`thermal_contact_heat` present". The pass
  and defect tests use an electrothermal pn diode built in the test file (an isothermal and an R_th sink, as U8); the
  three TVR-1 decks are NOT APPLICABLE (model off). The solver refuses energy balance, Schottky contacts and gate
  tunnelling with the model, so Σ V_c I_c covers every current-carrying contact.
- **C6 reading (owner ruling before code, 2026-11-04).** The solver forms each step's displacement current from the
  stored contact charges by that step's own formula (`src/solve/transient.cpp:598`: backward Euler after a start or
  breakpoint, variable-step BDF2 otherwise), and `error_ratio` is a max-norm of scaled ψ, densities and traps over
  rtol, not a charge: "the step's local error estimate summed" cannot be put in charge units from the results.
  Ruling: judged, per accepted step k and contact c, β h_k I_disp,k = Q_k − α₁ Q_{k−1} + α₂ Q_{k−2}, with the
  coefficients derived independently (BE: β = α₁ = 1, α₂ = 0; BDF2 with w = h_k / h_{k−1}: β = (1 + w)/(1 + 2w),
  α₁ = (1 + w)²/(1 + 2w), α₂ = w²/(1 + 2w)), the order and steps from the `TimePoint`s; tolerance 64ε times the sum of
  the operands' magnitudes. The cumulative sums are reported. Not judged: the trapezoidal ∫I dt against Q_N − Q_0
  per contact, with Σ_k error_ratio_k · rtol · \|ΔQ_k\| beside it, a diagnostic (as C9's Kramers–Kronig). Stated
  limit: the identity checks the bookkeeping (formula, order, restarts, scaling), not the charge physics, since the
  solver derives I_disp from the same stored Q.
- **C6 contacts (owner ruling before code, 2026-11-04).** The displacement identity on every contact; for each
  gate and meshed electrode also conduction current = 0 and terminal = displacement current, exactly. Predicate:
  a transient with a gate or electrode (else NOT APPLICABLE).
- **C6 test device (owner ruling before code, 2026-11-04).** `examples/moscap_cv.deck` loaded and run with
  `solve_transient` through the API (the deck is not edited): a gate waveform with a ramp and a step breakpoint
  (backward Euler restarts and variable-step BDF2), and a run with the backward-Euler integrator; defects injected
  on the result. The three TVR-1 decks are not transients: NOT APPLICABLE.
- **C7 tolerance (owner ruling before code, 2026-11-04).** §11.6's frozen bound as it stands: \|Y_ij − Y_ji\| ≤
  `resolution[k][i]` + `resolution[k][j]`, compared unnormalized (the same inequality as §11.6's ratio to max\|Y\|,
  without the division, so a zero matrix is judged). The resolution covers the AC solve's residual and the rounding
  of the current sums; it does not cover the operating point's distance from exact equilibrium (the Newton stopping
  error of C1 and C5) or the assembled matrix's rounding. No term is added for them before measurement; a miss is a
  finding for the owner. (The existing `tests/solve/small_signal_test.cpp:345` measured at most 0.06 of the bound.)
- **C7 models (owner ruling before code, 2026-11-04).** Every model is judged, no model predicate: Onsager
  reciprocity holds physically for each model at equilibrium. Discretely it is derived only for drift-diffusion with
  SG fluxes, recombination and traps (the linearized operator is symmetric in the quasi-Fermi variables); a miss under
  the electrothermal model, energy balance, Fermi–Dirac or gate tunnelling is a fail, recorded as a finding on the
  discretization (TVM2: nothing exempted in advance).
- **C7 devices (owner ruling before code, 2026-11-04).** With two contacts, zero column sums (C4) and zero row sums
  (gauge invariance) force Y_01 = Y_10, so a two-contact pass shows conservation, not reciprocity. The three TVR-1
  decks are run with a zero-bias AC point added through the API and a two-contact result carries a diagnostic saying
  so; the reciprocity evidence is a three-contact 2D MOS (gate, source, substrate; as `small_signal_test`'s `mos2d`)
  built in the test file: the pass at 0 V, and the same device biased judged with the predicate bypassed on a copy
  (the defect: a biased device is not reciprocal).
- **C7 shape (owner ruling before code, 2026-11-04).** `check_c7(SmallSignal, point, mesh, DeviceData)`: every
  frequency of the point (f = 0 included) as one part, the worst (frequency, i < j) named. NOT APPLICABLE when a bias
  is not exactly 0 or the thermal contacts are not at one temperature (C1's predicates); NOT EVALUATED when the
  resolution or admittance is missing or the point index is past the run's points. No magnetic field by
  construction (NiTCAD has none), as C1's generation.
- **C9 DC part (owner ruling before code, 2026-11-04).** The caller supplies four DC points about the AC point,
  contact j at V − 2h, V − h, V + h, V + 2h (every other bias equal). Conductance, every row i: Re Y_ij at f = 0
  against D4 = (8(I(+h) − I(−h)) − (I(+2h) − I(−2h))) / (12h); tolerance (8(r₋₁ + r₊₁) + r₋₂ + r₊₂) / (12h) (r the
  points' `terminal_current_resolution`) + \|D4 − D2\| (D2 = (I(+h) − I(−h)) / (2h), the stated truncation
  estimate) + 64ε Σ\|terms\| / (12h) + `resolution[0][j]`; f = 0 not among the frequencies: NOT EVALUATED.
  Capacitance, gate and electrode rows only (the DC results store no ohmic contact's charge; those rows NOT
  EVALUATED): Im Y_ij(f₁) / ω₁ at the smallest nonzero frequency against D4 of `gate_charge[i]`; tolerance \|D4 − D2\|
  of the charge + 64ε Σ\|terms\| / (12h) + `resolution[k₁][j]` / ω₁ + the ω² term \|C(f₂) − C(f₁)\| f₁² / (f₂² − f₁²)
  from the two lowest nonzero frequencies (fewer than two: NOT EVALUATED). The stencil not about the AC point
  (a bias off V + m h by more than 4ε\|V + m h\|, or another bias different): NOT EVALUATED. The DC points' Newton
  error in the charge is not in the bound (a miss is a finding).
- **C9 Kramers–Kronig (owner ruling before code, 2026-11-04; diagnostic only, never judged).** On C̃(ω) =
  (Y_ii(ω) − Y_ii(0)) / (jω) for each diagonal entry (causal in the e^{jωt} convention, finite at 0, → the
  geometric capacitance at ∞): Re C̃(ω) − C̃(∞) = −(2/π) P∫₀^∞ ω′ Im C̃(ω′) / (ω′² − ω²) dω′ (sign checked on
  1/(1 + jωτ)). Over the sampled band [a, b] the singular part is integrated in closed form and the regular part by
  the trapezoid rule in ln ω′; C̃(∞) is fixed at the highest interior sample (singly subtracted); the miss is
  reported at the other interior samples beside a tail estimate (2/π)(\|Im C̃(a)\| a² / (3ω²) + \|Im C̃(b)\| ω² /
  (3b²)) at ω and at the reference. Requires f = 0, ≥ 8 nonzero samples per decade and ≥ 4 decades; else the
  diagnostic says why and gives no value.
- **C10 limit continuity (owner ruling before code, 2026-11-04).** Generic: a caller-supplied x_k → 0 (≥ 3
  values, else NOT APPLICABLE "limit parameter swept (≥ 3 values)"), one model result per x_k and one limit-model
  result, a stated rate p. The difference e_k is one contact's \|I_k − I_lim\| (floor: the two resolutions) or
  max_i \|ψ_k,i − ψ_lim,i\| (floor: 64ε max \|ψ\|); any e_k not above its floor: NOT EVALUATED. The observed rates
  between successive x_k (metrics' `observed_orders_from_errors`) each within 0.2 of p. Demonstrations, rates
  derived here before measurement:
  (i) energy balance → drift-diffusion as τ_w → 0 with v_sat → ∞ (EC4's limit): with v_sat infinite Unit 52's
  mobility is μ₀ at any T_n, and the local balance gives T_n − T_L = (2/3)(q/k) τ_w μ₀ E² + O(τ_w²); the current's
  departure from drift-diffusion is the thermal diffusion ∇(n k T_n / q) less its lattice value, linear in
  T_n − T_L: p = 1 in τ_w.
  (ii) Fermi–Dirac → Boltzmann as N → 0: F_{1/2}(η) = e^η − e^{2η} / 2^{3/2} + …, so at neutrality n = N gives
  η_FD − η_B = N / (2^{3/2} N_c) + O(N²): the potential of a symmetric pn diode at 0 V moves by V_T N / (2^{3/2}
  N_c) on each side: p = 1 in N (the next term is N / N_c smaller).
  (iii) electrothermal → isothermal as κ → ∞ with every thermal contact an isothermal sink at the device
  temperature: the lattice heat equation is linear in κ for a fixed source, so T − T₀ = O(1/κ), and the current's
  change (mobility, n_ie, thermoelectric terms ∝ ∇T) is linear in T − T₀: p = 1 in 1/κ.
  MLDA → classical and Lombardi → bulk are listed for TVR-2/TVR-5 with their benchmarks.
- **C11 symmetry and invariance (owner ruling before code, 2026-11-04).** The caller supplies a node map and a
  contact map from result A (mesh, device) to result B. Predicates (else NOT APPLICABLE, "the map is a symmetry of
  the device", the first failing item named): the maps are bijections; doping, permittivity, region and insulator
  flags equal under the map, exactly; each edge of A maps to an edge of B with length within 64ε relative and
  coupling area within 64ε of the largest coupling at its ends; volumes within 64ε relative; contact kinds and node
  sets correspond; biases equal under the contact map, exactly. Tolerances: terminal currents \|I_A − I_B\| ≤
  res_A + res_B + 64ε max(\|I_A\|, \|I_B\|); gate charges 64ε max\|Q\| + the C3-style Newton term (Σ over the
  contact's edges C_f · 2 V_T, plus a lumped gate's C_ox A V_T, times u_A + u_B). Ruled before measurement
  (§11.6's 64ε relative alone is below two separate solves' Newton and resolution level). Devices: a 1D
  n+–n–n+ resistor mirrored (nodes reversed, contacts and biases swapped) and a 2D generated unstructured mesh
  rotated by 30° and translated, rebuilt by `make_dual_mesh`.
- **C8 passivity (owner rulings before code, 2026-11-04).** Per frequency k of a point, H = (Y + Yᴴ)/2 formed
  exactly Hermitian (H_ii = Re Y_ii, H_ji = conj H_ij), its eigenvalues by TC9's Jacobi; λ_min ≥ −tol with tol =
  √n ‖resolution[k]‖₂ (Weyl: an error E with \|E_ij\| ≤ r_j has ‖E‖₂ ≤ ‖E‖_F = √n ‖r‖₂) + 8nε‖H‖₂ (TC9's frozen
  accuracy) + the final off(A). Predicates C7's (every bias exactly 0, thermal contacts at one temperature);
  every frequency, f = 0 included; one part per point, the worst frequency named. More than 16 contacts (TC9's
  scope): NOT EVALUATED. Devices: C7's 3-contact mos2d copy and the three decks at equilibrium; defects at
  2 × tol and the bound's edge (0.99 / 1.01); with two contacts a diagnostic (with zero row and column sums
  λ = {0, 2 Re Y_00}, so C8 reduces to Re Y_00 ≥ 0). Y's null vector (1, …, 1) puts λ_min near 0 by conservation;
  the check is that it does not fall below −tol.
- **TVB-32 deck runs (owner rulings before code, 2026-11-04).** Scope: every check whose inputs the deck's own run
  supplies: C1–C5 and C12 at every operating point (the sweeps' points and the C–V run's DC points) and C1, C3 on
  the equilibrium solve; C7–C8 at every C–V point and at the zero-bias AC point added through the API (TVR1-Q5).
  C6, C9, C10 and C11 are listed as "input not produced by the deck" (their real-solve evidence is in their own
  tests). Acceptance: the expected outcome of each check at each point is fixed here, before any run; any other
  outcome fails the test and is a finding (TVM2). Frozen expectations (default models: no electrothermal, complete
  ionization, rates stored, no interface traps):
  C1, C2: pass where every bias is exactly 0 and the thermal contacts are at one temperature, else NOT APPLICABLE
  ("all biases 0"); C1 on the equilibrium: pass. C3: pass at every point and on the equilibrium (the Schottky and
  ohmic contact cells skipped by rule; the MOS-C gate a group). C4: pass. C5: NOT APPLICABLE ("electrothermal on").
  C12: pass. C7, C8: pass at a point whose biases are exactly 0, else NOT APPLICABLE ("all biases 0"); with two
  contacts each carries its two-contacts diagnostic.
- **Findings and departures of C5, C9, C10 and C11 (owner, "choose recommended", 2026-11-04).** None changes a
  tolerance (TVM2). C5: the heat rows' missing Newton term is accepted as a listed limitation; the hidden
  `[.finding]` case stays; option B (K₀ and the edge thermal conductances in `DeviceData`) is deferred, not
  refused. C9: the ω² term's blindness to an f₁-only error is accepted as a listed limitation of the ruled bound.
  C10: one limit result or one per value, accepted. C11: the 2D symmetric MOS mirror device is accepted in place
  of the 1D resistor; the 30° rotation alone stands in for the ruled rotation with translation; the generated
  mesh's corner tie (P2 refusal under rounding motions) is listed as a Unit 16 question, deferred.
- **TVB-01–03 (owner rulings before code, "choose recommended", 2026-11-05; TB-Q1–Q9).** Three misses were
  predicted before any run, from the code and its recorded measurements, and are frozen here as expected findings
  (no source outside the new targets changes, VR9; no tolerance moves, TVM2): (i) `base::hbar_J_s` is NIST's
  printed 10-digit truncation 1.054571817e-34 of the exact h/2π = 1.0545718176461565…e-34 (≈ 6.1e-10 relative;
  the header's "(exact)" is wrong); (ii) F₋½ = `fermi_half().derivative` is recorded at 3e-12 relative on
  (15, 40] (`tests/physics/fermi_dirac_test.cpp:79`) against TVB-02's 1e-12; (iii) `assemble::bernoulli` returns
  0 for x > 709.78 (std::expm1 overflows) where B(x) is a normal double up to x ≈ 715.
  TB-Q1: q, k, ε₀ and m₀ bit-equal to NIST's printed CODATA 2018 decimals; ħ = h/2π and k/q computed in
  double-double from the exact h, k, q and judged by TB-Q2; ħ misses (finding (i)). TB-Q2: the a-priori
  rounding bound γₙ = n u/(1 − n u), u = 2⁻⁵³: γ₃ for k in eV/K and for ħ (h's rounding, 2π's, the quotient's), γ₄ for
  V_T(T) at T ∈ {1, 4.2, 77, 300, 400, 1000} K. TB-Q3: source NIST's archived CODATA 2018 table, read before code
  (PA); SI-to-cm conversion by shifting the decimal exponent of the typed string (exact). TB-Q4: TVB-02
  reference from DLMF §25.12 (F_j(x) = Γ(j+1)⁻¹ ∫₀^∞ t^j / (e^{t−x} + 1) dt, PA), t = u², split at t = η when
  η > 0, tail cut with a written bound; primary tanh-sinh (TC4), second method adaptive Gauss–Kronrod (TC3),
  cross-checked against the polylogarithm series −Li_{j+1}(−e^η) (DLMF 25.12.16) for η ≤ −2; the reference
  counts only if its uncertainty is ≤ 1e-3 of the criterion. TB-Q5: Blakemore (paywalled) replaced by the exact
  anchors F_j(0) = (1 − 2⁻ʲ) ζ(j+1), ζ(3/2), ζ(1/2), ζ(5/2) typed from OEIS. TB-Q6: a hidden `[.finding]` case
  asserts ≤ 1e-12 everywhere; the visible case requires pass for F½, F₃/₂ and the inverse and that F₋½'s
  misses lie only at η > 15. TB-Q7: inverse reference η by Brent (TC5) in double-double on F_ref(η) = ν for each
  double ν; metric |Δη| / max(1, |η|) ≤ 1e-12; F₃/₂'s derivative and B′ recorded as diagnostics only. TB-Q8:
  321 η points (step 0.25 from −40 + 0.0123) plus −2 ± 1e-9, 40 − 1e-9 and 0; if Debug exceeds the remaining
  VR10 budget the Gauss–Kronrod cross-check runs at every 4th point (fixed now). TB-Q9: TVB-03 on a logarithmic
  grid of ±[1e-300, 745] including ±1e-2 and its neighbouring doubles; ≤ 4ε relative where B is a normal double,
  the subnormal range a diagnostic; reference x/expm1 in double-double, second method the series from exact
  Bernoulli numbers, scaled exponent (DD × 2^m) above x = 670; B on (709.78, ~715) misses (finding (iii)).
  Layout: `codata2018.hpp`, `fermi_dirac.{hpp,cpp}`, `bernoulli.hpp` in `nitcad_verify_reference`; tests in
  `tests/verify/benchmark/`, judged through `judge()` with filled contracts (VR3's first half).
- **TVB-02 finding (ii) widened (owner, "choose recommended", 2026-11-05, option A, after measurement).** The frozen
  "F₋½ passes for η ≤ 15" was not met: F₋½ misses from η = 13.76 (82 of 324 points, worst 2.0e-11 at η = 38.5), against
  the 3e-12 the header records. The expectation is not moved silently: the visible case now requires pass for F½, F₃/₂
  and the inverse and records F₋½'s outcomes; the hidden `[.finding]` case still asserts 1e-12 everywhere. No code
  change (VR9); the inaccurate accuracy statements of `physics/fermi_dirac.hpp` and the sampled existing gate
  (`tests/physics/fermi_dirac_test.cpp:79`) are listed for the owner.
- **Report and CLI (owner rulings before code, "choose recommended", 2026-11-05; RC-Q1–Q7).** RC-Q1: TVB-01/02/03/32
  move into `nitcad_verify_benchmark` as registered functions (id, tier, run → comparisons); the tests call them and
  assert the frozen outcomes; the TC1–TC9, metric and contract L0 tests stay in ctest and are named in the report as
  the reference's build-time evidence, not re-run. RC-Q2: a frozen known-findings table (benchmark, quantity,
  ruled region, DECISIONS reference) reports a `fail` inside a ruled region as `finding` (counted, exit 0, TVD-23);
  any other fail stays `fail`; a known finding that no longer reproduces is reported as pass with the note "known
  finding (x) no longer reproduces: remove the entry". RC-Q3: VR9 exception, narrow: `io::Content::verify_report = 7`
  with `save_verify_report` / `load_verify_report` (a new io header and source; existing io code and tests
  unchanged), plus a CSV summary beside it (`<report>.csv`). RC-Q4: the thread count in the report's own header, not
  in `results::Provenance`; VR7 compares two loaded reports with `started_utc` and `wall_s` cleared. RC-Q5: the three
  TVB-32 decks embedded at build time (a generated header), read by `deck::parse(text)`. RC-Q6: an unknown `--id` is
  usage (64) listing the known ids; `--tier full|nightly` run the same set as `fast` in TVR-1; no `--report` → the
  summary on stdout only. RC-Q7: VR8 as a ctest on the built exe: 0 for pass/finding/NA/NE, 64 usage, 1 for an
  unwritable report path and for an injected fail from a test-only registry entry `TVB-TEST-FAIL` compiled only into
  a test build, never into the shipped CLI.
- **VR10 (owner, "choose recommended", 2026-11-05, option A, after measurement).** The CLI-level VR7 test runs on
  TVB-01 and TVB-03 (`verify_cli.reproducible_TVB01`, `_TVB03`) instead of TVB-32, whose two deck runs (11.7 s in Debug)
  put the serial Debug verify run at 63.2 s against §11.9's 60 s. The budget is not moved. TVB-32's report was shown
  reproducible once (Release and Debug, 2026-11-05) and is no longer asserted at the CLI level; the in-process VR7
  case covers TVB-01 and TVB-03 as well.

**11.14 As-built (branch `verify/reference-numerics`, local).**

*TC1 — double-double (`nitcad_verify_reference`: `include/NiTCAD/verify/reference/double_double.hpp`,
`src/verify/reference/double_double.cpp`).* Sources read in full before code: Hida, Li and Bailey, "Library for
Double-Double and Quad-Double Arithmetic" (8 May 2008; Algorithms 3, 4, 7; §4.2 sqrt; §5.1 exp; §5.2 log) and
the QD library source (`include/qd/inline.h`, `include/qd/dd_inline.h`, `src/dd_real.cpp`, scibuilder/QD mirror)
for ieee_add (Briggs–Kahan), the FMA product and accurate_div, which the report does not print. Not read:
Dekker 1971 (paywalled), the ARITH-15 paper (its PDF fonts do not extract), Joldes–Muller–Popescu 2017 (HAL's
bot wall), so no published double-double error bound is cited; grade PB. Choices: sqrt by the report's Newton on
x⁻² (three steps; QD's source uses Karp's trick, the existing test helper one Newton step on x); exp with k = 256
(the report; QD's source uses 512); log by three Newton steps (the report; QD's source one) and, for \|a − 1\| <
1/8, 2 atanh((a − 1)/(a + 1)) by series so log is accurate relative to itself near its zero; π by Gauss's
formula 48 atan(1/18) + 32 atan(1/57) − 20 atan(1/239), ln 2 = 4 acoth 7 + 2 acoth 17 (derived in the source:
2 ln(4/3) + ln(9/8)), e = exp(1), √2 = sqrt(2). A first ln 2 formula written from memory (18 acoth 26 − 2 acoth
4801 + 8 acoth 8749) was wrong by 5e-11 and the constant test failed on it, as it should; the first ln 2
literal also had a dropped digit (a transcription error, caught by the same test and re-read against OEIS's
b-file).

Measured (Release, 10⁴ points each): constants π 0.30, e 0.016, ln 2 0.11, √2 0.24 u_dd (limit 1);
(√x)² = x 1.65 u_dd; (x / y) · y = x 0.86 u_dd; exp(log x) = x 0.34 u_dd per unit condition number;
exp(a + b) 0.81 per unit; added expm1(2x) identity 1.38 u_dd and expm1(log(1 + d)) = d 1.41 u_dd (limit 8 each);
two_sum, quick_two_sum, two_prod exact on all 10⁴ integer-scaled cases; added double-double + under cancellation
and × against the exact integer result, all within 8 u_dd. Added tests are additions only; no frozen tolerance
moved.

Mutations (each built and run, then reverted): two_sum without (b − v) — 7 tests fail; ieee_add without t.lo —
survived the identity tests (the term is ≈ 0.25 u_dd without cancellation), which is why the cancellation test
was added; it then fails; two_prod without the FMA error — at least 8 fail (output truncated at 8); × without a.lo·b.hi — 8 fail; ÷ keeping only q1
— 6 fail; ÷ without q3 — **survives** (q1 + q2 is already within ≈ 1 u_dd; detecting it needs a tolerance below
the frozen 8 u_dd, not set); exp squaring 7 times instead of 8 — 5 fail; sqrt with no Newton step — 2 fail; sqrt
with the Newton factor 1/4 — 2 fail (one Newton step instead of three is not a defect: from a double start one
step already converges quadratically, and it passes); log's atanh branch without the factor 2 — 2 fail; log with
one Newton step — 1 fails; π's series not alternating — 1 fails.

VR2 for this library: `tests/verify/boundary.cmake` runs as a build step before `nitcad_verify_reference`
compiles (and again as the ctest `verify_reference_boundary`); adding `#include "NiTCAD/physics/…"` to the
source, and separately `#include <Eigen/Core>` to the header, each failed the build at that step (shown once,
2026-11-04); `tests/verify/reference/header_boundary_test.cpp` checks Eigen's include path at compile time.
Ctest entries `verify.double_double`, `verify.header_boundary`, `verify_reference_boundary`: Release 0.16 s,
Debug 0.94 s.

*TC2 — cross-check with the test helper (`tests/verify/reference/cross_check_test.cpp`).* The test includes the
helper; the library does not (the boundary step would refuse it). Measured (Release, 10⁴ points each, TC1
grids): +, −, ×, ÷ 0 u_dd (bit-identical: the same published algorithms); sqrt 0.76; exp 0.42; expm1 1.65; log
0.25 absolute (after ruling 1); pow 0.93 per unit max(1, \|b\|, \|b ln a\|) (after ruling 2); π 0.16; ln 2 0
(limit 2 each). Mutations of the new library, each caught by TC2: sqrt Newton factor 1/4 (sqrt fails); exp
squaring 7 times (exp/expm1 and log/pow fail); log's atanh without the factor 2 (log/pow fails). A probe that
judged pow against repeated squaring was discarded: the squarings' own error doubles each step, so it was not
an adequate reference.

*TC3 — Gauss–Kronrod 7/15 (`include/NiTCAD/verify/reference/gauss_kronrod.hpp`,
`src/verify/reference/gauss_kronrod.cpp`).* Sources read before code: netlib `quadpack/dqk15.f` (data statements,
the error estimate) and `quadpack/dqage.f` (the adaptive driver), Piessens and de Doncker. The rule is computed, not
typed: Gauss nodes by Newton on P₇; Kronrod nodes as the roots of the Stieltjes polynomial E₈ = P₈ + Σ c_j P_j
(j = 0, 2, 4, 6) fixed by ∫ P₇ E₈ P_k = 0 (k = 1, 3, 5, 7; the integrals by a 16-point Gauss–Legendre rule, exact
for their degree 22), bracketed between the Gauss nodes and polished by Newton; Kronrod weights by exactness for
P₀ … P₁₄. The 33 printed digits are used only by the test. dqk15's estimate and dqage's loop are followed line by
line with ε = u_dd; deviations: the largest-error interval by linear search (not dqpsrt), the running error sum
in double-double (in double its drift, ≈ 2⁻⁵³ of the first errors, swamped a 1e-24 target and went negative —
the first run reported negative estimates), and the bad-integrand endpoint test in double-double.

Measured (Release): nodes ≤ 0.14 u_dd, Gauss weights ≤ 0.95, Kronrod weights ≤ 2.41 against the printed digits;
xᵏ on [−1, 1] within 1.64 u_dd for even k ≤ 22, 1.5e24 u_dd at k = 24, odd k exactly 0; ten adaptive integrals at
ε_rel = 1e-24 (eˣ, 1/(1+x), 1/(1+x²), √x, log(1+x) on [0, 1]; 1/(x + 1e-6); 1/(\|x − 0.3\| + 1e-4); T₄₀ and T₁₀₀ on
[−1, 1]; e^(−50x)): all status ok, true error ≤ estimate ≤ tolerance, 1 to 542 intervals (T₁₀₀), largest true error
3.4e-27 (the kink). Mutations, each caught: one Kronrod weight × (1 + 1e-20) (3 of 4 tests fail); a wrong
Stieltjes condition (k = 9 for 7; 3 fail); the error exponent 3 for 1.5, the raw estimate 1000× too small, and the
error sum rounded to double each step (the adaptive test fails on each). Ctest `verify.gauss_kronrod`: Release
< 0.1 s, Debug 0.17 s.

*TC4 — tanh–sinh (`include/NiTCAD/verify/reference/tanh_sinh.hpp`, `src/verify/reference/tanh_sinh.cpp`).*
Takahasi and Mori (1974) is paywalled and was not read; graded PB under TVR1-Q2 from two open reproductions read
in full before code, which agree on everything used: Bailey, Jeyabalan and Li, Experimental Mathematics 14 (2005)
317–329, §§3–5 (QUADTS), and Bailey, "Tanh-sinh high-precision quadrature" (19 Jan 2006), §§2–4. As they describe
it: x_j = tanh u₂, w_j = u₁ / cosh² u₂, u₁ = π/2 cosh hj, u₂ = π/2 sinh hj; level k uses h = 2⁻ᵏ from level 1,
later levels evaluate only the odd-indexed abscissae; 1 − x_j stored as 1 / (e^{u₂} cosh u₂) and [a, b] mapped
in double-double (their secondary-precision scaling); their heuristic estimate E = 10^d, d = min(0, max(d1²/d2,
2d1, d3, d4)), E = 1 for n ≤ 2, 0 when S_n = S_{n−1}, with ε₁ = u_dd and d not rounded (optional in both). Two
owner rulings after measurement (§11.13): ε₂ = ε₁³ (the note's n ≥ 2 family) and the 8 u_dd · ∫\|f\| rounding
floor. Infinite intervals are mapped by the caller, as BJL do: ∫₀^∞ e^{−x²} by x = 1/s − 1. The table (finest
level 10, 5 056 to 5 119 pairs: j·2⁻¹⁰ up to t ≈ 4.95) is computed once.

Measured (Release, ε_rel = 16 u_dd): ∫₀¹ x^{−1/2} = 2 0.49 u_dd (estimate 7.9e-31, true 4.9e-32); ∫₀¹ ln x = −1
0.89 u_dd (3.9e-31 / 4.4e-32); ∫₀¹ x^{−1/2} ln x = −4 1.47 u_dd (1.6e-30 / 2.9e-31); ∫₀^∞ e^{−x²} = √π/2 0.42 u_dd
(3.5e-31 / 1.8e-32); all status ok at level 4 (159 evaluations) except e^{−x²} at level 6 (633); added ∫₋₁¹ 1 and
∫₀¹ eˣ within 0.11 u_dd. Every estimate ≥ 2.6× the true error and ≤ ½ the tolerance. Added: ∫₀¹ x^{−0.9}, which
the table cannot resolve (tail ≈ 10 (r y)^{0.1} ≈ 2e-9), ends `limit_reached` at level 10 with estimate 5.0e-8 ≥
true error 2.2e-9 (d4 reports it); a loose tolerance still runs three levels. Mutations (each built, run,
reverted), all caught: w = u₁ / cosh u₂; y = 1 − tanh u₂ (cancels); u₂ from cosh t; ε₂ = ε₁² and ε₂ = 1e-40;
later levels re-adding the even points; h off by 2; the right point b + dx; the left point a + dx/2; the rounding
floor removed (the bound test fails, as measured before the ruling); 2d1 dropped. Two survived the first test set
and were caught after additions: d4 replaced by d3 (no frozen integrand makes d4 bind once n = 3; the x^{−0.9}
test catches it) and trusting levels 1–2 (no test had a loose tolerance; the loose-tolerance test catches it).
Ctest `verify.tanh_sinh`: Release 0.03 s, Debug 0.16 s (verify suite Debug 1.88 s).

*TC5 — Brent's zero (`include/NiTCAD/verify/reference/brent.hpp`, header only).* Source read before code:
Brent, *Algorithms for Minimization without Derivatives* (1973), chapter 4, in the author's scan
(maths-people.anu.edu.au/~brent/pd/rpb011i.pdf; image-only, rendered page by page): the ALGOL 60 procedure
*zero* (§6, pp. 58–59), the tolerance (2.9) δ = 2ε\|b\| + t, the error bound (2.18), the evaluation bound
(3.2)–(3.4) (k + 1)² − 2, and Table 4.1 (p. 55). PA. netlib `go/zeroin.f` was compared; it uses 2ε\|b\| + tol/2
and is not followed there. The translation is line by line in double-double with macheps = u_dd; additions:
refusal when f(a), f(b) have the same nonzero sign (zero "assumes" otherwise; zeroin.f checks), and a stop at
twice (3.4) so an excess is reported, not looped. Roots' digits: Wallis's number (OEIS A007493, 64 digits) and
the Dottie number (A003957, 68 digits), PA, compared exactly. cos is not in the library; the test uses a local
Taylor series, checked by cos² + sin² = 1 (1.06 u_dd) and cos x = 2cos²(x/2) − 1 (1.88 u_dd) on [−1, 1].

Measured (Release, t = u_dd): x³ − 2x − 5 on [2, 3] 0.0074 u_dd in 10 evaluations (bound 10 607); cos x − x on
[0, 1] 0.16 u_dd in 10 (bound 11 023); e^x − 1 on [−1, 2] \|ζ̂\| = 1.6e-13 u_dd in 13 (bound 11 447); an
interval without a sign change refused after 2 evaluations. Added tests: the smooth roots in ≤ ¼ of bisection's
k + 1 (10 of 103 and 105); Brent's own example (2.1), on which Dekker's procedure takes ≈ 256 steps, in 21
(bound 79); Table 4.1's x⁹ rows in 78 and 191 evaluations against the published 81 and 189 (IBM 360,
macheps = 16⁻¹³), asserted within ±5% — a band set after that measurement, to separate mutations the frozen
criteria do not; the x¹⁹ row is not used (x¹⁹ underflows in double below \|x\| ≈ 1e-16; Brent ran it with
zero2's extended exponent). Mutations (built, run, reverted): caught — the (b − a) sign in inverse quadratic
interpolation; no minimum step of tol; no swap keeping \|f(b)\| smallest; interpolation always accepted; m =
c − b; the bracket check removed; the bound as k² − 2; linear interpolation without its factor 2; the int/ext
branch inverted; e taken from the last step. Before the added tests only the last four of these (m, bracket,
bound, branch) were caught. **Survived:** forced bisection removed (both triggers, or only \|e\| < tol) — on
every test function the acceptance test p < \|s q\|/2 already rejects the crawling steps, and no tested function
reaches Brent's worst case; and the relative term 2ε\|b\| of the tolerance dropped — double-double resolves steps
far below u_dd·\|b\|, so termination and counts are unchanged here. Ctest `verify.brent`: Release and Debug
0.01 s (verify suite Debug 1.90 s).

*TC6 — Dormand–Prince 5(4) with dense output (`include/NiTCAD/verify/reference/dopri5.hpp`, header only).*
Sources read before code: Hairer and Wanner's DOPRI5 (Fortran, 25 April 1996, corrected 8 August 2005; the code
of HNW I, 2nd ed., §II.5), in the copy shipped with SciPy 1.10.0 (`scipy/integrate/dop/dopri5.f`; the authors'
page refused the download), all of DOPCOR, HINIT, CONTD5 and CDOPRI; and Boost.Odeint's `runge_kutta_dopri5.hpp`
(Ahnert and Mulansky), whose tableau and 4th-order weights agree with dopri5.f's E_i = b_i − b̂_i. Not read:
Dormand and Prince 1980 (ScienceDirect refused), Shampine 1986 (AMS bot wall), HNW's book. PB, as §11.13 rules.
DOPCOR is followed line by line in double-double with UROUND = u_dd (t and h double-double; the error norm and
step factors double); deviations: no stiffness detection, no SOLOUT call at the initial point, the last step lands
on tend exactly. A fixed-step driver (t_k = t0 + k h) serves the order tests.

Measured (Release): tableau against Boost 0 u_dd (c, a), e against b − b̂ 0.09 u_dd; row sums 0.63 u_dd; the 17
order-≤ 5 conditions of b 0.28 u_dd; b̂'s 8 order-≤ 4 conditions 0.15 u_dd and an order-5 residual 8.1e-4 (its
order is real); b(θ)'s 8 conditions at θ = 0.1 … 1 0.28 u_dd, b(1) = b exactly, end slopes 0.05 u_dd. Fixed-step
orders (N = 128 → 256 → 512 → 1024): y′ = y 4.990 / 4.995 / 4.998; oscillator on [0, 2π] 4.997 / 4.999 / 5.000;
y′ = y² 4.864 / 4.935 / 4.968. Continuous order (h = 2⁻³ … 2⁻⁷): y′ = y 3.86–3.98, oscillator 4.00, y′ = y² 4.11–4.01.
Proportionality slopes over 1e-14 … 1e-18: 0.998, 1.000, 0.991, every run `ok`, landing on tend exactly, nfcn =
6·nstep + 2. Global dense error at mid-steps (printed only): orders 4.87–5.11. Added tests: y′ = 5t⁴ exact to
0.20 u_dd and y′ = 6t⁵ not (7.5e-7), since no autonomous problem sees the stage times; a first step h0 = 1 rejected
4 times (not counted, as dopri5.f) with the final error 5.3e-15 against 5.2e-15 from HINIT's step; refusal of
rtol = atol = 0 and of t0 = tend; nmax = 3 → `too_many_steps`; a backward run; the dense output covering [t0, tend]
without gaps, once per accepted step.

Mutations (built, run, reverted; tests failing of 11): caught — a₅₂ = −25361/2187 (7); e₅ = −17254/339200 (5);
d₄'s numerator +1 (3); CONT(4) with +h k₇ (2); stage times t instead of t + c_i h (1, the t⁵ addition only);
no FSAL reuse in the adaptive loop (3); no last-step clip (2); a rejected step keeping h (1, `too_many_steps`); e₇
dropped (5); the fixed-step t_k rounded to double (1); θ₁ = θ (2); CONT(3) from k₇ (2). **Survived** (step-size
selection only; no test judges efficiency): the exponent 1/4 for 1/5 − 0.75β; the norm without its square root
(err² ≤ 1 ⇔ err ≤ 1); acceptance at err ≤ 10 (no accepted step here has err in (1, 10]); the scale without
\|y₁\|; Lund stabilization removed; HINIT returning 10⁻⁶. Ctest `verify.dopri5`: Release 0.03 s, Debug 0.30 s
(verify suite Debug 2.28 s).

*TC7 — Chebyshev collocation with domain splitting (`include/NiTCAD/verify/reference/chebyshev.hpp`,
`src/verify/reference/chebyshev.cpp`).* Sources read before code: Trefethen, *Spectral Methods in MATLAB*, chapter
6 (the author's sample-chapter PDF, rendered page by page: Theorem 7, (6.3)–(6.6), cheb.m, the printed cheb(1) …
cheb(5)), the online cheb.m and p13.m, the errata list; Boyd, *Chebyshev and Fourier Spectral Methods* (2nd ed.,
2000), §§22.1–22.4 (patching, (22.3)–(22.5)) and §7.7 (conditioning). PA. The differentiation matrix is cheb.m's
((6.5) off the diagonal, (6.6) on it), with two exact identities in place of cheb.m's arithmetic: x_j =
sin(π(N − 2j)/2N) and x_i − x_j = −2 sin((i + j)π/2N) sin((i − j)π/2N) (sin by an internal Taylor series on
[−π/2, π/2]); D² = D·D. The two-domain solve collocates −ε(2/L)²D²u = f at interior points, Dirichlet at the outer
ends and Boyd's patching with (22.5b) written for the flux ε_L u_L′ = ε_R u_R′; Gaussian elimination with partial
pivoting in double-double, after exact power-of-two row equilibration. The equilibration was added after the
first measurement (the unscaled solve left the u-continuity row's residual at 2e-28, rising to 2e-26 at N = 128,
because collocation rows carry entries of order N⁴); it is a solve-quality change, not a tolerance. U_ref: the
last two Chebyshev coefficients' magnitude per subdomain (heuristic, reported only).

Measured (Release): cheb(1) … cheb(5) within 3.6e-5 of the printed 4-decimal values; nodes 1.22 u_dd from the
test's own cos; the (6.6) diagonal 0.23 of its bound; D·xᵏ 0.12 of its bound; Chebyshev coefficients of T_k 0.66
u_dd; Clenshaw on T₅ 1.25 u_dd; a 6 × 6 integer system 0.58 u_dd. The ε-jump problem, N = 4 … 64: interpolant
error 5.93e-2, 1.24e-5, 2.49e-15, 4.74e-29, 3.35e-28 (floors 5.1e-29 … 1.3e-26); u(0) the same; left flux
5.91e-2 … 2.26e-29, 3.72e-28, right 2.29e-28, 3.63e-28 (floors 1.5e-27 … 3.9e-25); at N = 64 the u mismatch 0 and
the flux mismatch 7.3e-28. p13.m's problem split at 0: 6.08e-2, 3.94e-5, 8.10e-14, 4.07e-30, 3.34e-28 (floors
1.1e-28 … 2.7e-26). Invalid problems (N < 2, c outside (a, b), no f) refused.

Mutations (built, run, reverted; tests failing of 6): caught — the diagonal left 0 (4); c_j/c_i (4); x_i − x_j's
sign (3); nodes from 2N + 1 (4); the coefficient end-halving dropped (3); Clenshaw's last step with 2ξ (3); no
pivoting (1); the flux row's sign (2); the continuity row's sign (2); the left scale not squared (2); ε_left in
the right collocation rows (1); no equilibration (1, the u-continuity mismatch); the sin series stopped at 2⁻⁶⁰ (4);
the left flux with the right scale (1); the right piece evaluated with the left coefficients (2). **Survived:** the
tail estimate from one coefficient instead of two (the estimate is reported, not judged). Ctest
`verify.chebyshev`: Release 0.08 s, Debug 0.68 s (verify suite Debug 3.18 s).

*TC8 — complex-step differentiation in double (`include/NiTCAD/verify/reference/complex_step.hpp`, header only).*
Squire and Trapp (1998) and Martins, Sturdza and Alonso (ACM TOMS 2003) refused automated download (WVU repository
and ResearchGate 403, Academia.edu walled) and were not read. Graded PB under TVR1-Q2 on the same authors' two
precursors, read in full (aero-comlab.stanford.edu/Papers, rendered page by page): Martins, Kroo and Alonso, AIAA
2000-0689 ((11)–(14); the Lyness–Moler function (15) at x = 1.5; the step range to ≈ 1e-307; (19)–(20) and Table
5), and Martins, Sturdza and Alonso, AIAA 2001-0921 ((1)–(9): the arcsin and sin pitfalls, the linearised
definitions). Neither prints the derivative's digits at 1.5, so the reference is the analytic derivative in
double-double (owner ruling). `complex_step(f, x, h = 1e-200)` returns Im f(x + ih)/h, Re f(x + ih), and the
Richardson truncation estimate 1.25·\|d(2h) − d(h)\|/3 (two evaluations); it refuses non-finite x and h that is not
a positive normal double or whose 2h overflows, and reports `non_finite` when f does. `cs::abs`, `max`, `min`,
`sign`, `dim` are Table 5's new definitions (branches on real parts; abs takes +z at 0; sign's magnitude by (20),
where the Table writes \|x₁\|). The test's double-double sin and cos reduce by kπ/2 with the library's π and sum
Taylor series on \|r\| ≤ π/4.

Measured (Release, h = 1e-200, ε = 2⁻⁵²): the composite at 1.5 0.40ε (f′ = 4.05343…), its value 0.81ε (addition);
eˣ on [−200, 200] 1.58ε (1 000 judged); sin on [−2π, 2π] 0.98ε (994 judged; 6 near cos x = 0 printed, 0.24ε);
x^{3/2} as z·√z on [1e-50, 1e50] 1.24ε (1 000 judged); the composite on [0, 1.5] 2.42ε (1 000 judged). Printed
only: pow(z, 1.5) 113ε at x = 2e-50 (the exp–log condition number, as predicted before measuring). Richardson at h
= 1e-3: estimate/true = 1.25 (3 printed digits) for all four; at 1e-200 the estimate is 0. Observed order of
the truncation (composite, h = 2⁻⁴ … 2⁻⁶) 1.996, 1.999 (band 2 ± 0.1, an addition). Test-local sin, cos: sin² + cos²
− 1 1.0 u_dd, sin 2x − 2 sin x cos x 1.13 u_dd (limit 4). Table 5: abs −1 / +1 / +1 (right derivative at 0) at x =
−2, 2, 0; the modulus std::abs gives 0; max, min (ties to z₁), sign and dim give their branch derivatives. Refusals:
h = 0, −1e-200, NaN, ∞, 1e-310, DBL_MAX and x = NaN, ∞; h = DBL_MIN accepted. Martins' lower limit shown: eˣ at x
= −700 returns status ok with derivative 0 (e⁻⁷⁰⁰·1e-200 underflows inside f; the primitive cannot detect it, as
the header says). The 2001 pitfalls reproduced: arcsin by −i log(iz + √(1 − z²)) and sin by (e^{iz} − e^{−iz})/2i
both return a derivative far from exact; MSVC's std::asin agrees with 1/√0.75 at 0.5 to the 6 printed digits (not judged).

Mutations (built, run, restored; test cases failing of 8): caught — derivative Im/(2h) (5); Richardson /4 for /3,
the second step at h, d(2h) divided by h, F_s = 1 (1 each, the bound); abs with x ≤ 0 taking −z, max and min ties
to z₂, sign without (20), dim without its branch (1 each); subnormal h accepted, 2h overflow accepted, no
non-finite check (1 each); value = 0 (1, after the value check was added). F_s = 2 first **survived**: the upper
check used the library's own constant, so it scaled with the mutation; it now uses the ruled 1.25 and fails it.
**Survived:** the value taken from the 2h evaluation — at h = 1e-200 Re f(x + 2ih) = Re f(x + ih) exactly, an
equivalent mutant at the default step. Ctest `verify.complex_step`: Release 0.02 s, Debug 0.05 s (verify suite
Debug 3.15 s).

*TC9 — small Hermitian eigenproblems by cyclic Jacobi (`include/NiTCAD/verify/reference/jacobi_eigen.hpp`, header
only).* Source read before code: Forsythe and Henrici, Trans. AMS 94 (1960), in a scan (pages.mtu.edu; rendered
page by page): §1 (the cyclic orders (III_r), (III_c), condition (IV), the unitary form (5)–(6), Theorems 1–4), §2.1
(Lemma 1, (10)–(14)), §2.5 ((41), (42)) and §3.1 (Theorems 5 and 6, (43)–(50)); not read: the proofs of §§2.2–2.4
and §3.2 (singular values, not used). Golub and Van Loan (book) and Fettis and Caslin 1967 (AMS 403) not read;
Demmel and Veselić's LAWN 15 was downloaded but not read (not cited by TC9). PA on Forsythe and Henrici. The
rotation is FH (6) with α = δ = 0, β = θ, γ = −θ (a_pq = r e^{iθ}), φ from (44) as the smaller root t = tan φ, so
\|φ\| ≤ π/4 ((IV) with J = [−π/4, π/4], Theorem 3); cyclic by rows; a zero a_pq skipped; a′_pq set to 0 and the
diagonal to its real part after each rotation; stop at off(A) ≤ ε‖A‖_F; U_ref = the final off(A) (Weyl; rotation
rounding not included). A first version guarded τ² overflow with t = 1/(2τ) for \|τ\| > 1e150; no test could reach
it and without it t = 0 (error below rounding), so it was removed. The Hilbert coefficients are in
`tests/verify/reference/hilbert_charpoly_oeis.hpp`, generated from the b-file by awk (no digits typed); λ_min by
the first sign change upward from 1e-20 on a 200-per-decade grid, λ_max downward from 3.2, each polished by Brent
(t = 1e-45) in double-double.

Measured (Release; units nε‖A‖₂, limit 8): L·H_n, n = 3 … 10 (3–6 sweeps): λ_min ≤ 0.012, λ_max ≤ 0.89, residual ≤
0.92, orthogonality ≤ 2.24 nε (limit 8); λ_min(H_n) = 9.67023e-5 (n = 4) … 1.09315e-13 (n = 10), relative errors
4.3e-15 … 1.3e-4 (printed only: far better than the absolute bound, as Demmel and Veselić's relative-accuracy
result for definite matrices would predict, not judged here). D ± zzᴴ, n = 8 and 16 (5 sweeps): eigenvalues ≤ 1.00,
residual ≤ 1.02, orthogonality ≤ 3.36; I + zzᴴ (n = 7, the multiple eigenvalue) 0.08; M D Mᴴ, n = 6 and 16 (5, 6
sweeps) ≤ 0.71; FH (42) 0.04 in one sweep; a diagonal matrix exact, sorted, 0 sweeps; refusals of n = 0, n = 17, a
size mismatch, a non-Hermitian pair, an imaginary diagonal (1e-300), NaN; the 16 × 16 ones matrix accepted.

Mutations (built, run, restored; test cases failing of 8): caught — τ's sign flipped (5); the larger root, \|φ\| > π/4
(4); the phase factors swapped (3); V not rotated in p (5); V with the wrong phase (3); A's columns with the wrong
phase (3); stop at 1e-8 (3); no sort (5); no Hermitian check (1); n = 17 accepted (1); a one-sweep cap (3); zero
a_pq not skipped (3, NaN). **Survived:** a′_pq not set to 0, and the diagonal's imaginary residue kept — both
residues are at rounding level and the next sweep removes them (eigenvalues are read from the real part). With
TC9, VR1's L0 set (TC1–TC9, each with a recorded mutation) is complete. Ctest `verify.jacobi_eigen`: Release 0.02
s, Debug 0.12 s (verify suite Debug 3.17 s).

*Metrics (§11.5; `include/NiTCAD/verify/reference/metrics.hpp`, `src/verify/reference/metrics.cpp`, in
`nitcad_verify_reference`).* Source read before code: NPARC "Examining Spatial (Grid) Convergence" with VERIFY's
`prD.do`/`prD.out` (PA for the formulas printed in text and the worked example; the page's equation images carry
no content of their own). In double, finest grid first: relative error with floor and log error with floors
(`not_evaluated` at or below a floor); L∞, dual-volume L2 and dual-face H1 (the value returned is √Σ c_e Δe²; §4.2
writes the sum, its square), summed in double-double with exact power-of-two scaling; orders from errors (any
ratio) and from three solutions (constant r); Richardson; GCI (F_s 1.25 / 3.0, relative to the finer solution, a
fraction); `grid_convergence` over N ≥ 3 grids, refusing (`not_evaluated`, with the reason) non-constant r (1e-12
relative), a zero difference, non-monotone solutions over the whole set, any p̂_k ≤ 0 and a zero solution; the
asymptotic check per the ruling (≥ 4 grids: spread of successive p̂ ≤ 0.1 and ρ_k = GCI_{k+2,k+3}(p̂_k) /
(r^{p̂_k} GCI_{k+1,k+2}(p̂_k)) within 10 % of 1; ρ_k = r^{p̂_{k+1} − p̂_k} f_{k+1}/f_{k+2}, so it is not the
three-grid identity; the implementer's choice of pairs, open to the owner); ε_int summed exactly in double-double
before one rounding; U = Σ\|U_k\| + \|ε_int\|, interacting when \|ε_int\| > 0.1 U.

Measured (Release): NPARC p̂ = 1.7861696, f_{h=0} = 0.9713003 (increment 0.0008003), GCI₁₂ = 0.1030826 %, GCI₂₃ =
0.3562493 %, each within half a unit of the printed digit; the page's ratio recomputed from its rounded
intermediates 1.0020189 (prints 1.002019); the study's three-grid ratio equals f₁/f₂ = 1.0020237 to 4ε; the
asymptotic check `not_evaluated` ("needs 4 grids"). Synthetic orders p ∈ {1, 1.5, 2, 4} (errors, r ∈ {2, 1.5, √2,
1.1} and a non-constant set) and p ∈ {1, 1.5, 2, 3} (solutions) recovered to 1e-12; Richardson recovers f* to
1e-12. Corrections to my own test inputs, before any tolerance was considered (none moved): a power-law set
f = 1 + 1e-5 h² cancelled to 1.3e-12 in p̂ and was replaced by the exact dyadic f = 1 + 2⁻²⁰h²; an expected
ρ deviation named the wrong pair (1 − f₃/f₄ for 1 − f₂/f₃); a "ratio out of range" input with ρ = 0.917
(miscomputed as 0.84) was replaced by C = 2⁻⁵ (ρ = 0.75).

Mutations (31, built, run, restored; test cases failing of 10): first run 25 caught, 5 survived, 1 did not build
(an unused loop variable under /WX; rewritten as an always-false condition, then caught). Survivors, each a test
gap closed by a new assertion and then caught: the stand-alone non-monotone refusal removed (the NaN order still
refused it through p̂ ≤ 0; the test now checks the reason); the order-spread test and the ratio test each ignored
(no case had one failing alone; one of each added); ε_int accumulated in double (no cancelling input; 1 + 2⁻⁶⁰ − 1
added); L∞ without abs (no negative maximum; added). Final: 31 of 31 caught. Ctest `verify.metrics`: Release
< 0.01 s, Debug 0.01 s (verify suite Debug 2.97 s). VR4 met.

*Contract and outcomes (§11.4; `include/NiTCAD/verify/reference/contract.hpp`, `src/verify/reference/contract.cpp`,
in `nitcad_verify_reference`, namespace `NiTCAD::verify`).* Plain data: `ReferenceClass` (E–X), `Independence`
(I0–I3), `Provenance` (PA–PD), `Outcome` (pass, finding, fail, NOT APPLICABLE, NOT EVALUATED), `Predicate`,
`Quantity`, `ReferenceInfo` (the §3.3 fields, with `second_method` for E-with-PD's two methods and E-with-I1's
cross-check) and `Comparison`; `validate` for each, returning the problems by name. `judge(info, predicates, met,
criterion)` applies, in order: the contract (invalid → NOT EVALUATED, problems named); predicate validity; every
listed applicability predicate supplied (else NOT EVALUATED, the name given); every supplied predicate holding (else
NOT APPLICABLE, "name = value (threshold t)", whether or not the criterion was met); then I0 → finding "regression
pin, not evidence"; PC → finding "provenance insufficient"; else pass or fail. `missing_input` builds a NOT
EVALUATED with the input named. Numbers in reasons are the shortest round-trip form (`std::to_chars`; a first
`%g` loop printed 700 as "7e+02", caught by the test). The I0 rule is applied before PC, so an I0 PC reference
reports the I0 reason. Choice: a predicate the caller supplies beyond the contract's list still refuses when false
(the consistency checks add their own); PC remains a valid grade to hold (the E rule makes it invalid for E only).

Tests (7 cases, 137 assertions): a filled contract passes or fails on the criterion; 17 single-field corruptions
each give exactly one named problem and NOT EVALUATED; the E/A rules satisfied by what they ask for (PD with a
second method, I1 with a cross-check, A with an error, PB for M); applicability refusals; **PC cannot pass or fail
(VR3)**; I0 likewise, applicability still first; quantity, missing-input and outcome validation. Mutations (20,
built, run, restored): 20 of 20 caught (one first failed to build under /WX, an unused variable, and was rewritten
as an always-false condition). Ctest `verify.contract`: Release 0.01 s, Debug 0.02 s (verify suite Debug 3.96 s,
13 tests). VR3's second half is met here; its first half (every TVR-1 reference carries a filled contract) is met
when TVB-01–03 and the TVB-32 checks are written with their contracts.

*Consistency, DC checks C1, C2, C4, C12 (§11.6; `include/NiTCAD/verify/consistency/{device_data,checks}.hpp`,
`src/verify/consistency/checks.cpp` as `nitcad_verify_consistency`; `verify/benchmark/device_data` as
`nitcad_verify_benchmark`).* **Deviation from §11.2's table, an implementer's choice under the owner's "choose
recommended" (2026-11-04), open for review:** `device::Device`'s header includes `physics` headers (its regions
hold `physics::Semiconductor`), so a consistency layer that included `device` would reach `physics`
transitively. The layer reads a plain `DeviceData` (doping, permittivity, insulator flags, contacts, whether
interfaces carry traps or fixed charge) that the benchmark layer copies from the device; its allowed includes are
verify/consistency, verify/reference, results, mesh and base. VR2 for this layer: `tests/verify/boundary.cmake`
before it compiles (direct includes) and `tests/verify/consistency/header_boundary_test.cpp` (transitive: one type
of physics, assemble, solve, linalg and device is first declared as a variable of the same name, so any path
that reaches that layer's header fails to compile; the Eigen include-path check is not made there because the
test links the deck). Shown once each, restored: `#include "NiTCAD/device/device.hpp"` in `checks.hpp` stopped the
build at the boundary step; `physics/semiconductor.hpp` placed before the consistency headers in the test gave
C2061/C4430 at `semiconductor.hpp:452`.

Each check fills a contract (`check_info`: class D, I2, PD, the §11.6 row as citation and derivation) and is
judged through `judge`: applicability first (NOT APPLICABLE with the predicate and value, e.g. "all biases 0 =
0.3 (threshold 0)"), then parts (each with its worst residual, tolerance and location); any failing part fails the
check; else any NOT EVALUATED part makes the check NOT EVALUATED (a check is not passed on the parts it could
evaluate alone); else pass. C1: zero-bias and one-thermal-temperature predicates; terminal currents against their
own resolutions; each edge's electron and hole current against the smallest nonzero contact resolution (the
reading of "edges ≤ that resolution", a choice); the quasi-Fermi part as ruled above. C2: each stored mechanism's
q Σ\|R\|V in double-double against the smallest nonzero resolution; NOT EVALUATED with no rate stored. C4: \|Σ I\| ≤
Σ resolution, and for the weighted current when present; NOT EVALUATED without currents. C12: n, p > 0 on
semiconductor nodes, every stored temperature (lattice, carrier, thermal contacts, device) > 0, exact; trap
occupancies are not in the results, so that part is NOT EVALUATED when an interface has traps.

Measured on `examples/pn_diode_1d.deck` (anode 0 V and 0.3 V set through the API; Release): C1 at 0 V pass —
terminal currents 2.2e-6 of their resolution (2.80e-6 A/cm²), edge currents 2.6e-3 (worst electron edge 123),
quasi-Fermi 0.51 of the ruled bound at node 112; the equilibrium solve's C1 3.5e-4 (node 191); C1 and C2 at 0.3 V
NOT APPLICABLE with the bias named; C2 at 0 V pass (SRH 1.4e-19, Auger 9.4e-23, radiative 0 of the resolution);
C4 pass at both points (1.1e-6, 1.3e-4); C12 pass. Injected defects, each failing its check at the named
location: a terminal current at twice its resolution (contact anode), one hole edge current at 1.01 × the bound
(edge 150), F_n shifted 1e-12 eV (node 200) and, on the equilibrium result, F_p shifted 1e-13 eV (node 50); a
recombination of 2 × the resolution at one node, and a ±2 × pair that cancels in a signed sum; a current
corrupted by 3 Σ resolution, a NaN current; a negative, NaN or zero density (n and p), a zero temperature. NOT
EVALUATED shown for missing bands, a missing Newton record, no stored rates, missing currents and interface traps.

Mutations of `checks.cpp` (24, built, run, restored; one malformed pattern dropped before the run): first run 20
caught, 4 survived, each a test gap closed and then caught: C2 without \|·\| and without the volume (the injected
defect failed either way; a ±pair and a below-threshold rate added), C12 with n ≥ 0 (no zero electron density;
added), a NaN residual passing (no NaN residual reached a comparison; a NaN current added). Final: 24 of 24
caught. Ctest `verify_consistency.dc_checks`: Release 0.02 s, Debug 0.37 s (verify suites Debug 4.01 s, 16 tests).

*Consistency, C3 dual-cell Gauss law (§11.6 and the two C3 rulings of §11.13; `check_c3` for a `BiasPoint` and an
`EquilibriumResult` in `src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/gauss_test.cpp`).*
Restated from the discretization as read (`src/assemble/scaled_device.cpp`, `terms.cpp`, `gate.hpp`,
`interface_edges.hpp`): the Poisson variable is `potential_V` / V_T (band shifts enter the bands, not the flux); an
edge's coupling is ε₀ ε_f A/d with ε_f the harmonic mean of its ends' permittivities; a semiconductor row's charge is
V (p − n + N_D − N_A); a lumped gate adds G(ψ_G − ψ) + q Q_f A per node and stores Σ G(ψ_G − ψ) as `gate_charge`; a
charged interface's edge is two half-edges whose fluxes sum to the sheet charge q σ_f A + trapped (stored per
interface, without σ_f). `DeviceData` gained per-node regions, per-interface region pair, `charged`, σ_f and total
trap density (levels plus D_it × band width), per-gate face areas, C_ox = ε₀ ε_ox / t_ox and Q_f, and the run's
`incomplete_ionization` flag (`benchmark::device_data(device, models)`). Groups by disjoint sets: each contact's
nodes, each charged edge's ends, each charged interface's edges. Predicate "converged state" (else NOT APPLICABLE);
NOT EVALUATED with incomplete ionization, missing fields, no Newton record, or gate or trap charges not stored
(then the groups and the global part both). Choices: a charged edge adds 2 ε_max (A/d)(\|ψ_a\| + \|ψ_b\|) to the
rounding scale (its half-edge conductances; ψ_I is not stored), a group's stored charges add their magnitudes, a
lumped gate adds 2 C_ox A \|ψ\| to the scale and C_ox A V_T to the Newton term, a trapped interface 2 q N_t A
(occupancy moves by at most 2u for an update u); the p-poly gate's electrothermal shift is not in the Newton term
(a listed gap: no TVR-1 deck has it).

Measured (Release; worst residual / tolerance): pn diode 0 V and 0.3 V pass (worst cell 0.0019, 0.0059); pn equilibrium (0.0084)
pass; MOS capacitor at gate 0, −0.5, −1, 0.25 V pass (cells ≤ 0.0067, the gate group ≤ 0.0071); Schottky diode at
0, −0.5, 0.2 V pass (≤ 0.0083), the anode's charge from its cells' flux −3.9e-8 C/cm² at 0 V; a 1D silicon/oxide
stack written in the test (electrode, interface σ_f 5e11, a trap level and band) at gate 0, ±0.05 V and
equilibrium pass (interface group ≤ 1.0e-3); the same stack without the interface (harmonic mean in use), a moscap
with Q_f 3e11 under its gate, and a 2D oxide device (11 × 9, many charged edges, from `build_test`) at equilibrium
pass (≤ 7.9e-4). Injected defects, each failing at the named place: ψ + 1e-9 V (node 200), n and p moved 1e-6
(node 150), gate charge dropped, an extra or dropped Q_f, interface trap charge dropped (1D and 2D), σ_f dropped,
electrode charge dropped, ψ + 1e-9 V on a 2D interface node; a contact cell's corruption is not judged (skipped).
A test restates the ruled tolerance independently and holds the reported one to it (1e-9 relative) for a cell
(with and without u_last), the lumped gate and the trapped interface.

**Finding (not chased, for the owner):** the 1D silicon/oxide stack's DC sweep with a 0.25 V step on the electrode
stops at the second point: "Newton iteration 3: pivot ratio below min_pivot_ratio: the matrix is singular to working
precision"; ±0.05 V steps converge. The test uses the converging points; the deck is in the test file.

Mutations of the C3 code (25, built, run, restored): first run 17 caught, 6 survived, 2 malformed by the runner (a
`||` in the pattern split its `|`-separated fields; rewritten without it, then caught). Survivors, closed by the
ruled-tolerance test and then caught: the rounding scale without the charge term, one end's face scale, the
interior-edge skip (it changes the Newton term), one end's face Newton term, the trap Newton term, the gate oxide
Newton term. Final 25 of 25 caught. A first version of that test kept a pointer into a temporary result; it passed
in Release and crashed in Debug, and was fixed before the final mutation run. Ctest `verify_consistency.gauss`:
Release 0.07 s, Debug 1.62 s (verify suites Debug 5.59 s, 17 tests).

*Consistency, C5 steady energy balance (§11.6 and the three C5 rulings of §11.13; `check_c5` for a `BiasPoint` in
`src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/energy_test.cpp`).* \|Σ_c V_c I_c − Σ_k H_k\|,
summed in double-double, against E_max Σ_c R_c + 64ε (Σ_c \|V_c I_c\| + Σ_k \|H_k\|) as ruled. Predicates
"electrothermal on" (`DeviceData::electrothermal`, the run's model) and "`thermal_contact_heat` present" (else NOT
APPLICABLE); NOT EVALUATED when biases, currents, resolutions, one heat per thermal contact, the bands, the lattice
temperature or the thermopower exponents are missing. `DeviceData` gained the model flag and per-node thermopower
exponents r_n, r_p (0 on insulator nodes), copied by `benchmark::device_data`. The per-node heat is a listed gap
in the contract's assumptions, not a part. Contract: class D, I2, PD, equations G9 and G15.

Measured (Release) on an electrothermal 1D pn diode written in the test (1 µm, n 1e17 | p 1e17 by an erfc step,
41 nodes, isothermal sink at x_min, R_th 1e-3 K cm²/W at x_max, anode on the p side at x_max; default Newton
tolerance): 0, 0.3, 0.5, 0.7, 0.8 V pass at 2.7e-4, 9.7e-5, 1.6e-4, 4.3e-4, 7.0e-4 of the bound (0.8 V: 1.31 kW/cm²,
residual 4.2e-11 W/cm², bound 5.9e-8). `examples/pn_diode_1d.deck` at 0.3 V is NOT APPLICABLE ("electrothermal
on"), as are the other two TVR-1 decks (model off). Injected defects, each failing: a sink's heat ±2 × the bound,
the R_th sink's heat dropped, the anode current moved by 2 × bound / V, the heats' sign flipped, a NaN heat, and
with the resolutions zeroed (the rounding term alone) a 1e-12 relative imbalance. A test restates the ruled
tolerance independently and holds the reported one to it, with each term moved alone (r_n, r_p with \|r + 5/2\|,
T_max, an insulator node's band skipped, a contact voltage above every band energy); another holds the benchmark
copy to silicon's r = −1/2 on every node.

**Finding (for the owner): the ruled C5 bound omits the heat rows' Newton stopping error.** On the same bar fully
compensated (N_D = N_A = 1e17, a resistor, R_c ≈ 5e-15 A/cm², so E_max R ≈ 1e-14 W/cm²) with the anode at x_min
next to the isothermal sink, at the default tolerance (1e-8) the 0.8 V point misses: residual 2.89e-14 against
1.61e-14 W/cm² (1.8 ×), last Newton update 1.30e-9; 0.5 and 0.7 V reach 0.51 and 0.40 of the bound (updates 2.1e-9,
1.3e-9). At tolerance 1e-12 the same points end with updates near 5e-15 and sit at 4.9e-4, 3.9e-4, 6.7e-4. By U1's
identity (Σ K₀ f_T = −Σ V I at any state) the imbalance is the heat rows' residual, which E_max R covers only where
R is large (U8's devices). With the anode at x_max the resistor passes (≤ 1.5e-3). The tolerance is unchanged
(TVM2); the hidden case `[.finding]` "Finding C5: the compensated resistor at the default Newton tolerance" reports
both tolerances and fails at the miss. Option B of the question (a Newton term from u_last through the heat rows)
would cover it; that needs K₀ and the edge thermal conductances in `DeviceData`.

Mutations (26, built, run, restored): 24 of `check_c5` (predicates, each NOT EVALUATED input, signs of power and
heat, each magnitude, resolution and E_max term, r with and without 5/2 and \|·\|, T_max, the insulator skip, the
rounding term, a signed residual, a doubled bound) and 2 of the benchmark copy (model flag, r_n); first run 26 of
26 caught. Ctest `verify_consistency.energy`: Release 0.02 s, Debug 0.36 s (verify suites Debug 6.08 s, 18 tests;
Release 18/18).

*Consistency, C6 transient charge bookkeeping (§11.6 and the three C6 rulings of §11.13; `check_c6` for a
`results::Transient` in `src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/transient_test.cpp`).*
Per accepted step k ≥ 1 and contact c, β h_k I_disp − (Q_k − α₁ Q_{k−1} + α₂ Q_{k−2}) in double-double against
64ε (\|Q_k\| + \|α₁ Q_{k−1}\| + \|α₂ Q_{k−2}\| + β h_k \|I_disp\|), the coefficients derived in the check (BE; BDF2 from
the quadratic through the last three times, w = h_k / h_{k−1} from the stored steps). Gates and meshed electrodes
without gate tunnelling: \|I_cond\| + \|I_term − I_disp\| = 0 at every point, exactly (`ContactData::tunnelling` added,
copied by `benchmark::device_data`). NOT APPLICABLE without a gate or electrode; NOT EVALUATED with fewer than two
points, charges or currents not one per contact, an order other than 1 or 2 (or 2 without two earlier points), or a
step not positive. `CheckResult` gained `diagnostics` (name, value, scale, where, note; never judged): per contact
the cumulative identity and the trapezoidal ∫I_disp dt − (Q_N − Q_0) beside Σ error_ratio · rtol · \|ΔQ_k\| (rtol
from the run record's `transient.rtol`; NaN with a note when absent; steps without an estimate counted). Contract:
class D, I2, PD, equations G1 and G15, the bookkeeping limit in its assumptions.

Measured (Release) on `examples/moscap_cv.deck` through `solve_transient` (gate 0 → −0.5 V over 100 ps, held,
jump to 0.1 V at 200 ps; body 0; t_end 400 ps; rtol 1e-2 in the test to keep the Debug cost down, VR10): BDF2 run
160 points (6 BE, 153 BDF2, 2 BE restarts after BDF2), identity worst 0.26 of its bound; backward-Euler run 251
points, worst 0.0067; the gate part exactly 0 on both. At the default rtol (1e-3) the runs had 280 and 679 points
(worst 0.43 and 0.0072) and C6 took 35 s in Debug. Diagnostics: the cumulative identity is at rounding (≤ 1e-21
against \|Σ ΔQ\| ≈ 6e-9 C/cm²). **The trapezoidal integral is off by about a quarter of the gate's charge change**
(−1.6e-9 against Q_N − Q_0 ≈ 6e-9 C/cm² with BDF2; −5.2e-10 with BE; the summed error scale 7e-10 to 9e-10), the
body nearly the opposite. The likely cause, not chased: the 0.6 V jump at 200 ps. The step after it (backward
Euler) carries the whole jump of Q as I_disp = ΔQ / h, which a trapezoid halves. That is the reason the judged part
uses the step's own formula; the diagnostic is reported as measured. A pn-diode transient
(`examples/pn_diode_1d.deck`, anode ramp) is NOT APPLICABLE ("transient with a gate or electrode").

Injected defects, each failing its part: the gate's I_disp at the first BDF2 step moved by twice the bound (a
1e-9 relative change there sits below the rounding of Q, since that step's ΔQ is small against Q: the first form of
this defect passed and was resized against the bound), the body's at the first BE step 1e-9, a charge 1e-9 (also at
the last point), a BDF2 step relabelled BE, the previous step 1e-6 longer, a NaN charge; a gate conduction current
of 1e-300, a gate terminal current one ulp off its displacement current; a tunnelling gate's conduction current is
not judged. Meshed electrode: the gate relabelled an electrode is judged the same way.

Mutations (25): 24 of `check_c6`, first run 20 caught; survivors closed by tests and then caught: the electrode
kind in the predicate (relabelled electrode), the β h \|I\| term in the bound (a two-point case whose residual 100ε
lies between 64ε\|Q\| and the full bound), skipping the last point (a defect there), rtol ignored in the diagnostic
scale (scale proportional to rtol). Final run on the final tests (rtol 1e-2): 24 of 24 caught. The benchmark copy of the tunnelling flag survives: equivalent on every device
here (no tunnelling gate in TVR-1), listed. Ctest `verify_consistency.transient`: Release 0.68 s, Debug 17.1 s
(verify suites Debug 22.8 s, 19 tests; Release 19/19).

*Consistency, C7 reciprocity at equilibrium (§11.6 and the four C7 rulings of §11.13; `check_c7` for one point of
a `results::SmallSignal` in `src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/reciprocity_test.cpp`).*
At every frequency k of the point and each pair i > j, \|Y_ij − Y_ji\| (the real and imaginary differences exact in
double-double, then the modulus) against `resolution[k][i]` + `resolution[k][j]`; one part, the worst frequency and
pair named ("f = 1e+08 Hz, Y(source, gate)"). NOT APPLICABLE at a bias not exactly 0 or thermal contacts at more
than one temperature (C1's predicates, from the point's `dc`); NOT EVALUATED for a point index past the run's
points, the run's contact count not the device's, admittance or resolution not one per frequency or not sized to
the contacts. With two contacts the result carries the diagnostic "two contacts" (the worst residual and bound,
and the note that zero column and row sums force Y_01 = Y_10). Contract: class D, I2, PD, equation G14, the
assumptions as ruled (no magnetic field; every model judged; the resolution's two omissions; the two-contact limit).
`checks.hpp` now includes `results/small_signal.hpp` (inside the layer's boundary; VR2's check passes).

Measured (Release; frequencies 0, 1e2, 1e5, 1e8, 1e11, 1e13 Hz; worst \|Y_ij − Y_ji\| / (res_i + res_j)):
the three-contact 2D MOS (`small_signal_test`'s `mos2d` copied into the test: a meshed metal gate, a source on an n+
pocket, the substrate; fixed charge, two trap levels, a trap band and surface recombination at the interface) at 0 V
**0.245** at f = 0, Y(source, gate) (`tests/solve/small_signal_test.cpp:345` measures the same 0.245 and 102 today;
its comment's "0.06" and "114" are stale from Unit 22, not edited here, listed for the owner); the largest off-diagonal \|Y\| at 1e13 Hz is about 5e12 times the largest resolution, so the check is not
vacuous. The same device at (−0.4, 0, 0.1) V is NOT APPLICABLE; judged with the biases relabelled 0 on a copy it
fails at **102 ×** the bound (100 Hz, Y(source, gate)): a biased device is not reciprocal. The three TVR-1 decks
with a zero-bias AC point added through the API (their own options; the decks not edited): `pn_diode_1d` 0.044,
`moscap_cv` 0.235 (1e8 Hz, Y(body, gate)), `schottky_diode_1d` 0.029, each with the two-contact diagnostic.

Injected defects, each failing: Y(source, gate) moved by i·2(res_0 + res_2) at 1e8 Hz (named); Y(gate, substrate)
by 2(res_0 + res_1), real, at f = 0 (named); a NaN admittance. A symmetric change of 5(res_0 + res_2) to both
entries passes. The bound itself: on a symmetrized copy, an asymmetry of 0.99 (res_i + res_j) at the pair and
frequency with the most nearly equal resolutions (1e13 Hz, (source, substrate), the smaller 0.498 of the sum) passes
and 1.01 fails.

Mutations (20, built, run, restored): the point-index guard; each predicate; each NOT EVALUATED input (contact
count, admittance count, resolution count, inner sizes); the bound as res_i alone, res_j alone, doubled; adjacent
pairs skipped; Y_ij against itself; the real or the imaginary part dropped (first written as an unused variable,
which `/WX` refused to build, then rewritten as 0 × the part); the names swapped; the frequency label; the
diagnostic dropped and always added; the first and the last frequency skipped. 20 of 20 caught on the first built
run. Ctest `verify_consistency.reciprocity`: Release 0.05 s, Debug 0.99 s (verify suites Debug 19.0 s, 20 tests;
Release 20/20).
Remaining in §11.6: C8–C11 and the TVB-32 runs on the three decks (C7's deck runs are in its test).

*Consistency, C9 time–frequency consistency at low frequency (§11.6 and the two C9 rulings of §11.13; `check_c9`
with a `DcStencil` in `src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/low_frequency_test.cpp`).*
Two judged parts. Conductance, every row i of the stencil's column j: \|Re Y_ij(0) − D4\| against (8(r₋₁ + r₊₁) +
r₋₂ + r₊₂)/(12h) + \|D4 − D2\| + 64ε Σ\|terms\|/(12h) + `resolution[0][j]` (the differences of currents in
double-double). Capacitance, gate and electrode rows: \|Im Y_ij(f₁)/ω₁ − D4(Q)\| against \|D4 − D2\| + 64ε
Σ\|terms\|/(12h) + `resolution[k₁][j]`/ω₁ + \|C(f₂) − C(f₁)\| f₁²/(f₂² − f₁²), f₁ < f₂ the two lowest nonzero
frequencies (sorted, any order given). With ohmic rows present the diagnostic "ohmic rows' capacitance" says they
are not evaluated; with no gate the part is NOT EVALUATED and so is the check (a check is not passed on part of
itself). NOT EVALUATED besides: a point index past the run's, admittance or resolution not sized, the stencil's
contact, step or points invalid or not at V + m h (4ε) with the other biases equal, f = 0 absent (conductance),
fewer than two nonzero frequencies (capacitance), currents or gate charges not stored. Kramers–Kronig as ruled,
per diagonal entry, as a diagnostic (value: the worst miss; scale: the tail estimate at f and f_ref); not estimated
(with the reason) without f = 0, with fewer than three nonzero samples, under 4 decades or under 8 per decade (the
widest step between neighbours). (Y − Y(0))/(jω) is written as (b − ja)/ω: MSVC's complex division warns of a
divide by zero under `/WX`.) Contract: class D, I2, PD, equation G14, no applicability predicate (a stencil not
about the point is a missing input, not a predicate: the contract requires every listed predicate supplied).

Measured (Release). `examples/moscap_cv.deck` at gate −0.5 V (depletion), stencil on the gate, h = 2 mV,
frequencies 0 and 100 Hz – 1 THz at 8 per decade: conductance worst 2.4e-6 of its bound (Y(body, gate));
capacitance C(100 Hz) = 1.2565982309863e-7 F/cm² against D4 1.2565982309375e-7 (3.9e-11 relative), bound 4.4e-12
(3.5e-5 of C), dominated by resolution/ω₁. From 1 Hz the bound was 3.7e-10 (3e-3 of C): the resolution's rounding
part does not fall with ω, so f₁ is chosen where resolution/ω₁ is small (recorded in the test). Kramers–Kronig on
the MOS-C: worst miss 8.3e-10 F/cm² (6.6e-3 of C) at 5.6 GHz against a tail estimate 4.6e-10 (f_ref 750 GHz).
`examples/pn_diode_1d.deck` at 0.5 V, stencil on the anode, h = 1 mV: conductance worst residual 1.9e-6 against
6.1e-3 S/cm² (G = 4.80 S/cm²; the bound 1.3e-3 of G, dominated by the DC resolution term ~1/h, balanced near
h = 1 mV by \|D4 − D2\| ~ h²: 1.6e-2 at 0.3 mV, measured); capacitance NOT EVALUATED (no gate), so C9 NOT
EVALUATED. The AC resolution is 1.2e-5 of the diode's conductance bound; on the MOS-C gate row (no DC current,
resolution 0) it is the whole bound. Synthetic Debye admittance G₀ + jω(C∞ + ΔC/(1 + jωτ)), 1e-4/τ – 1e4/τ at 8 per
decade: Kramers–Kronig miss 1.9e-3 ΔC; a real bump of 0.1 ΔC at one sample is seen at 0.10 ΔC.

**Limit of the ruled bound (for the owner):** the ω² term is estimated from C(f₁) itself, so an error in the f₁
sample alone raises the bound by f₁²/(f₂² − f₁²) = 1.29 (8 per decade) times the error and is never caught: a
change of 3 × the bound at f₁ alone passes (tested and recorded as passing). An error common to the low frequencies
(every Im Y scaled) is caught. Defects (each failing): Re Y(body, gate) at 0 Hz by 2 × the body row's bound (named);
Im Y(gate, gate) scaled to move C(f₁) by 3 × the bound; a DC gate charge at V + h and a DC current at V − h moved by
6h × the bound (D4 moves 8/12h, \|D4 − D2\| at most 2/12h of a change: residual ≥ 3 against a bound ≤ 2). The bounds
at their edge: 0.99 passes and 1.01 fails for the diode's anode conductance (the bound restated in the test), the
MOS-C gate row's conductance (the AC resolution alone) and the gate capacitance (Im Y scaled).

Mutations (29, built, run, restored): D4's weight and denominator, D2 for D4, each term of both bounds (the
truncation, the resolution weights, the AC resolution, the ω² term, resolution/ω₁ without ω₁), imaginary for real,
the transposed entry, C from Re Y, C at 2ω₁, the f = 0 test, the stencil offsets, its slack (1e10 ε), the null-point
and contact guards, the sort of the frequencies, every row as a gate, the two-frequency, charge and ohmic-diagnostic
guards; Kramers–Kronig's sign, its d ln ω, the decade and per-decade rules, Y(0) not subtracted, the conversion's
sign, the closed-form singular part dropped. First run: 20 caught, 6 survived, 3 refused to build (an unused
variable under `/WX`; rewritten as 0 × the term). The survivors were closed by tests and then caught: the bound
edges (truncation, resolution weights, AC resolution), a step 1e-9 relative off (slack), a contact past the device's
with every point at the AC biases (the guard), G₀ in the synthetic admittance (Y(0) subtraction). **One survives:**
the singular part dropped. Near-equivalent on the Debye test: P∫₀^∞ dω′/(ω′² − ω²) = 0, so the closed form is small
away from the band edges, where this Im C̃ is small; a diagnostic, listed. 28 of 29 caught. Ctest
`verify_consistency.low_frequency`: Release 0.09 s, Debug 2.75 s (verify suites Debug 26.7 s, 21 tests; Release
21/21).

*Consistency, C10 limit continuity (§11.6 and the C10 ruling of §11.13; `check_c10` with a `LimitSeries` in
`src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/limits_test.cpp`).* The values sorted
ascending; per value the difference (a contact's current, or the potential's L∞, in double-double) against its
floor; the observed rates between neighbours from `reference::observed_orders_from_errors`, each \|p − p_stated\|
against 0.2, the worst pair named ("tau_w [s] = 1e-15 .. 1e-14 (rate …)"); the differences and floors reported as
diagnostics. NOT APPLICABLE under 3 values; NOT EVALUATED (each guard naming itself) for model or limit counts, a
result missing, a value not finite and positive or repeated, the stated rate not finite and positive, the contact
out of range, currents or potentials not stored, a difference not above its floor, or a rate the metrics library
refuses. **Departure from the ruling's wording, for the owner:** the series takes one limit result *or one per
value*. The Fermi–Dirac → Boltzmann demonstration needs it, since its limit model's device moves with the doping
(the ruling's "one limit-model result" fits demonstrations (i) and (iii) only). Contract: class D, I2, PD,
equations G1–G3, applicability "limit parameter swept (>= 3 values)".

Measured (Release; rates derived in §11.13 before measuring, all stated 1):
(i) energy balance → drift-diffusion, 1 µm diode 1e17/2e17 (101 nodes) at 0.6 V, v_sat = 1e30 for both carriers,
τ_w = 1e-14, 1e-15, 1e-16 s against drift-diffusion with the default silicon: \|ΔI\| 1.42e-3, 1.68e-4, 1.71e-5 A/cm²
(floors 5.4e-7); rates 0.929 (1e-15 .. 1e-14) and 0.992. Pass; the coarser pair's 0.07 is the next order showing.
(ii) Fermi–Dirac → Boltzmann, symmetric diode N/N, 10 µm (201 nodes), 0 V, N = 1e15, 1e16, 1e17: max \|Δψ\| 3.6361e-7,
3.6361e-6, 3.6359e-5 V (floors 4e-15 to 6e-15); rates 1.000000 and 0.999978. The coefficient implies a density of
states of 2.514e19 cm⁻³ in V_T N/(2^{3/2} N_DOS) (reported, not judged; not cross-checked against the material's N_c,
N_v).
(iii) electrothermal → isothermal, 1 µm diode 1e17/1e17 (41 nodes) with isothermal sinks at 300 K at both ends,
0.7 V, κ scaled by 1, 10, 100 (x = 1/scale) against the isothermal device: \|ΔI\| 1.1508e-4, 1.1508e-5, 1.1508e-6
A/cm² (floors 6e-8); rates 0.999996 and 1.000000.

Defects, each failing: stated rates 0.75 and 1.25; the middle value's difference doubled (rates 1.30 and 0.70); the
largest value's alone doubled (only the last pair moves, to 1.23, named 1e-15 .. 1e-14); the FD study's first point's
largest node difference doubled. A difference of 1.8 × one resolution (below the two together) is NOT EVALUATED.

Mutations (18, built, run, restored): the predicate's 3, each guard, the sort, the per-value limit, the floor's
second resolution, the floor test, the potential's maximum (last node only), the 0.2 widened and narrowed, a signed
rate difference, only the first rate judged, the null guard. First run 14 caught; four guards survived because the
metrics library or a later size check refuses the same input with another reason (repeated values, x ≤ 0, the
contact) or nothing reaches them (the rate); closed by tests that assert each guard's own reason and a zero rate.
18 of 18 caught. Ctest `verify_consistency.limits`: Release 0.05 s, Debug 1.25 s (verify suites Debug 28.8 s, 22
tests; Release 22/22).

*Consistency, C11 symmetry and invariance (§11.6 and the C11 ruling of §11.13; `check_c11` with a `SymmetryMap` in
`src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/symmetry_test.cpp`).* The predicate "the map is
a symmetry of the device" is evaluated as ruled, in this order, the first failing item named in the NOT APPLICABLE
reason: node and contact counts; both device data consistent with their meshes; the node map, then the contact map,
a bijection; per node doping, permittivity, region and insulator flag exactly, the volume to 64ε relative; edge
counts; each of A's edges has an image (by its mapped ends), its length to 64ε relative and its coupling to 64ε of
the largest coupling at its four ends (A's and B's); contact kinds and node sets; biases exactly. Two judged parts:
terminal currents (\|I_A − I_B\| in double-double against res_A + res_B + 64ε max\|I\|, the worst named
"contact drain -> source"); gate charges, only with a gate or meshed electrode, against 64ε max\|Q\| + (Σ over the
contact's edges ε₀ ε_h A/d · 2V_T, plus a lumped gate's C_ox A V_T) (u_A + u_B), NOT EVALUATED without a Newton
record or stored charges. Contract: class D, I2, PD, applicability "the map is a symmetry of the device".

**Departures from the ruled devices, for the owner.** (1) The mirror device is a 2D symmetric MOS structure, not
the ruled 1D n+–n–n+ resistor: oxide on [−5 nm, 0), p-silicon 1e17 to 120 nm, n+ pockets 1e19 under the source
(y < 0.75 µm) and the drain (y > 2.25 µm) to 30 nm, width 3 µm (7 lines), a meshed metal gate (4.4 eV) on all of
x_min, a fixed interface charge 3e11 cm⁻²; contacts gate, substrate, source, drain; A at (1.0, 0, 0, 0.1) V, B at
(1.0, 0, 0.1, 0) V, contact map {0, 1, 3, 2}, node map y → W − y by coordinates. It exercises what the resistor
cannot: a gate charge, four contacts and a contact map that is not the identity. (2) **The ruled rotation by 30°
with a translation cannot be run on a generated mesh: `make_dual_mesh` refuses it** ("boundary facet is not Gabriel
(P2): node 53 is strictly inside its diametral sphere"), and so does any translation alone. The generated mesh of
[0, 2e-4] × [0, 1e-4] cm (`max_size` L/12) has a square cell at the corner (L, H): nodes 19 (the corner), 20 and 21
(along the two sides) and 53 (the diagonal) are exactly cocircular, so 53 lies exactly on the diametral circles of
the boundary facets (19, 20) and (19, 21) (measured: (a − p)·(b − p) = 0 in double-double; the predicate returns 0,
a tie, which P2 accepts). Moving the nodes by any motion that rounds the coordinates breaks the tie, and the
translation (3.7e-4, −1.3e-4) and the rotation with it break it inward. Exact motions keep it: the quarter turns
(−y, x), (−x, −y), (y, −x) are accepted; a dyadic translation is not exact on the generated interior nodes. The
30° rotation about the origin alone is accepted, so the test runs that: an observation about this mesh's rounding,
not a property of the motion. This is not a predicate defect (the predicates are exact; P2 is right about the
moved coordinates); it is that generated meshes are Gabriel only by a tie at domain corners, so rigid motion of a
generated mesh is not robust under `make_dual_mesh`. A Unit 16 question for the owner (tie-free corners in the
generator, or a moved mesh repaired rather than refused); no Unit 16 code was changed.

Measured (Release). Mirror: currents worst residual 1.6e-14 A/cm against 6.8e-12 (drain → source; the drain current
is above 1e3 × the bound), gate charges bit-identical (residual 0 against 2.2e-26 C/cm, whose Newton term is
restated in the test and matched to 1e-12 relative). Rotation, the generated pn diode (acceptors 1e17 for x < L/2,
donors 1e17 above, ohmic contacts on x = 0 and x = L) forward at 0.4 V: I = 2.4995e-8 A/cm, residual 1.4e-17 against
1.07e-12. Both pass. NOT APPLICABLE, each with its item named: one node's donors (×(1 + 1e-12)), permittivity,
region, insulator flag; a bias off by 1e-12 V; a contact map {0, 1, 2, 2}; two nodes mapped to one; the identity
node map (contact nodes do not correspond); a contact kind; B's mesh rebuilt by `Mesh::from_parts` with one volume
(×(1 + 1e-12)), one edge length (×(1 + 1e-13), inside `from_parts`' 1e-12) or one coupling (×(1 + 1e-12)) changed;
the unchanged mesh rebuilt passes. Defects: B's source current moved by twice its pair's bound fails, named; the gate
charge by twice the part's bound fails; no Newton record leaves the charges and the check NOT EVALUATED.

Mutations (24, built, run, restored): each predicate item disabled (doping, permittivity, region and insulator
together and each alone, volume, node and contact bijection, length, coupling, contact kind, contact nodes, bias);
the current bound ×4 on the resolutions, B's current at A's contact index; the charge part's Newton term dropped,
V_T for 2V_T, only edges with both ends on the contact, the arithmetic for the harmonic mean, u_B dropped, 64ε ×4;
no part for gates; the reason's item dropped; the Newton-record guard. First run (22): 15 caught, 4 survived
(region, length, contact kind, harmonic mean), 3 refused to build (unused variables under `/WX`; rewritten as
0 × the operands). The survivors were closed by the from_parts, region, insulator and kind tests, and the Newton
term by its restatement. **One survives:** the harmonic mean, equivalent here: every edge at the meshed gate joins
two oxide nodes, so ε_h is ε_ox either way. 23 of 24 caught. Ctest `verify_consistency.symmetry`: Release 0.07 s, Debug 1.83 s (verify suites at -j 6: Debug 18.7 s, 23
tests; Release 23/23).

Remaining in §11.6: C8 (passivity) and the TVB-32 runs on the three decks.

*Consistency, C8 passivity at equilibrium (§11.6 and the C8 ruling of §11.13; `check_c8` in
`src/verify/consistency/checks.cpp`; tests `tests/verify/consistency/passivity_test.cpp`).* Per frequency: H formed
exactly Hermitian (H_ii = Re Y_ii; H_ij = (Y_ij + conj Y_ji)/2, H_ji its conjugate); `reference::jacobi_eigen`;
residual max(0, −λ_min) against √n ‖resolution[k]‖₂ + 8nε max\|λ\| + off(A); one part, the worst frequency named
("f = 1e+08 Hz"). A non-finite H or resolution term fails (it does not reach the eigensolver, which would refuse it
as invalid input). NOT EVALUATED: a point index past the run's, the run's contacts not the device's, no contacts,
more than 16, admittance or resolution not one per frequency or not sized to the contacts, the eigensolver not
converged (unreachable on the tests' inputs: n ≤ 16, finite, Hermitian by construction). Diagnostics: the smallest
eigenvalue over the frequencies with its tolerance and frequency; with two contacts the note that the pass shows
Re Y_00 ≥ 0. Contract: class D, I2, PD, equation G14, C7's two applicability predicates. Forming H rounds each entry
once (≤ ε\|H_ij\|, so ≤ √n ε ‖H‖₂ in norm), inside the ruled 8nε‖H‖₂.

Measured (Release). C7's mos2d at 0 V, frequencies 0, 1e2, 1e5, 1e8, 1e11, 1e13 Hz: worst residual 7.1e-14 against
1.24e-12 (0.057) at 0 Hz. **Vacuous at 0 and 100 Hz (recorded):** jωC is symmetric and drops out of (Y + Yᴴ)/2, and
the equilibrium conduction (gate current 0, the source junction at zero bias) is below the resolution: λ_max 8.8e-13
against the resolution term 1.24e-12 at 0 Hz (0.71), 0.89 at 100 Hz. From 1e5 Hz λ_max is 1.7e4, 4.9e6, 2.0e9 and
5.3e10 times the resolution term; the test pins both bands. Decks at a zero-bias AC point (frequencies as C7's):
`pn_diode_1d` 1.0e-7 against 8.7e-6 at 1e5 Hz (0.012); `moscap_cv` 1.13e-7 against 5.46e-7 at 1e8 Hz (0.21);
`schottky_diode_1d` 1.6e-9 against 8.3e-8 at 1e8 Hz (0.020); each with the two-contacts diagnostic. On the decks the
smallest eigenvalue is negative at those levels (column and row sums zero only to the resolution). Defects, each
failing at the named frequency: Re Y(gate, gate) lowered by 6 × the tolerance at 1e8 Hz, at 0 Hz and at 1e13 Hz (along
the null vector u = (1, 1, 1)/√3, uᴴHu falls by δ/3); a NaN entry. A real skew change of 1e3 × that (Y_20 + d,
Y_02 − d) leaves H unchanged and passes. The bound at its edge on a synthetic Y = aP − c uuᵀ (eigenvalues −c, a, a;
a = 10 × the resolution term): c = 0.99 × the bound passes, 1.01 fails; with the resolution zeroed, Y = aP (λ_min 0
to rounding) passes.

Mutations (24, built, run, restored, each confined to `check_c8`'s lines): √n dropped; the 8nε term dropped; off(A)
dropped; the resolution's 1-norm for its 2-norm; λ_max for λ_min; the residual 0; H without the conjugate; H from the
lower triangle alone; the diagonal kept complex; each predicate forced true; the 16-contact, no-contact, eigensolver
status and finiteness guards; the two-contacts diagnostic at three; the smallest-eigenvalue tracking reversed; the
first and the last frequency skipped; the point, contact-count, resolution-size and frequency-count guards; the
diagnostic's tolerance as the resolution term alone. First run: 21 caught, 3 survived; the diagnostic's tolerance
was closed by a test (the defect's diagnostic equals its part) and then caught. **Two survive, argued:** off(A)
dropped (Jacobi stops at off ≤ ε‖H‖_F ≤ √n ε‖H‖₂, below the 8nε‖H‖₂ beside it, so no input at the edge tells them
apart); the eigensolver's status guard (non-convergence in 50 sweeps is not constructible for n ≤ 16, and every other
refusal is guarded before the call). 22 of 24 caught. Ctest `verify_consistency.passivity`: Release 0.06 s, Debug
1.26 s (verify suites at -j 6: Debug 18.8 s, 24 tests; Release 24/24).

Remaining in §11.6: the TVB-32 runs on the three decks. C1–C12 are built.

*TVB-32, the consistency suite on the three decks (§11.6, TVR1-Q5 and the TVB-32 ruling of §11.13; tests
`tests/verify/consistency/decks_test.cpp`).* Each deck read and built as `nitcad run` does, its outputs cleared
(nothing written), run by `deck::run`; the device data from `benchmark::device_data` with the plan's models. Every
operating point gets C1–C5 and C12, the equilibrium solve C1 and C3, every C–V point and the added zero-bias AC point
(`solve_small_signal` with the deck's own frequency, 1 MHz) C7 and C8. Each result is held to the frozen
expectation; a failed or NOT EVALUATED part also fails the test. C6, C9, C10 and C11: not run, input not produced by
the deck.

Measured (Release), every outcome as frozen, no finding:

| Deck | Results | C1 | C2 | C3 | C4 | C5 | C12 | C7 | C8 |
|---|---|---|---|---|---|---|---|---|---|
| `pn_diode_1d` (sweep 0–0.7 V, 15 points) | 15 + equilibrium | 2 pass (0 V, equilibrium), 14 NA; worst 0.51 (0 V, quasi-Fermi, node 112) | 1 pass, 14 NA; 1.4e-19 | 16 pass; 0.0085 | 15 pass; 4.9e-4 | 15 NA | 15 pass | — | — |
| `schottky_diode_1d` (sweep −1–0.4 V, 29 points) | 29 + equilibrium | 2 pass, 28 NA; 0.0096 | 1 pass, 28 NA; 7.4e-18 | 30 pass; 0.0097 | 29 pass; 3.1e-4 | 29 NA | 29 pass | — | — |
| `moscap_cv` (gate −2–2 V, 17 points, 1 MHz) + the added AC point | 18 + equilibrium | 3 pass, 18 NA; 0.013 | 2 pass, 18 NA; 1.5e-15 | 21 pass; 0.0082 | 20 pass; 1.1e-13 | 20 NA | 20 pass | 2 pass, 18 NA; 0.17 | 2 pass, 18 NA; 0.0056 |

The sweeps' grid values hit 0 V exactly (−1 + 20 × 0.05 and −2 + 16 × 0.125), so the zero-bias points are judged,
not NOT APPLICABLE by rounding. Every zero-bias C7 and C8 carries its two-contacts diagnostic. Ctest
`verify_consistency.decks`: Release 0.16 s, Debug 5.9 s. Verify suites: Release 25/25, Debug 25/25; serial Debug
30.7 s, within §11.9's 60 s budget.

§11.6 complete: C1–C12 built and TVB-32 run. Open for the owner: the listed limitations of §11.13's
"choose recommended" ruling; C3's 0.25 V electrode-step singular pivot finding; the stale comment at
`tests/solve/small_signal_test.cpp:345`. Next in §11: TVB-01–03 (§11.8) and the report/CLI (§11.7).

*TVB-01, constants and thermal voltage (§11.8, TB-Q1–Q3; `include/NiTCAD/verify/reference/codata2018.hpp`,
`tests/verify/benchmark/constants_test.cpp`, new executable `nitcad_verify_benchmark_test`, ctests
`verify_benchmark.*`).* Source read before code: NIST's `allascii_2018.txt` (header "2018 CODATA adjustment"); the
seven lines typed as printed, digit groups and "..." kept. The reference converts the printed digits itself:
`shifted_decimal` moves the exponent in the string (F m⁻¹ → F cm⁻¹ by −2), `to_double` is `std::from_chars`, and
`to_dd` forms M·10^e with 10^n exact as a double-double (10²² · 10^(n−22) by `two_prod`, n ≤ 44), one rounding.
Contract: E, I3, PA, G1 and G4.

Measured (Release), every outcome as frozen. q, k, ε₀ (F/cm), m₀: bit-equal to the printed decimals. k/q [eV/K]
1.58e-16 relative (γ₃ = 3.33e-16). V_T at 1, 4.2, 77, 300, 400, 1000 K: 1.58e-16, 1.83e-16, 1.20e-16, 1.14e-16,
1.14e-16, 1.14e-16 (γ₄ = 4.44e-16). **Finding (i), confirmed:** `base::hbar_J_s` is 6.13e-10 relative below
h/2π (bit-equal to NIST's printed truncation 1.054571817e-34, a diagnostic); outcome fail, as frozen; the hidden
`[.finding]` case asserts γ₃ and fails. Not changed (VR9): the fix (the literal 1.0545718176461565e-34, or
h / 2π in base) and its effect on band-to-band, gate tunnelling and quantum confinement are for the owner.
Reference self-checks: exponent shift and parsing, 10^n exact for 23 ≤ n ≤ 44, `to_dd` rounding to `from_chars`'s
double, the truncated printed ħ and k/q below their exact values by less than one unit of their last digit.
Mutations (9, built, run, restored): a digit of q, k, ε₀, m₀ or h; 10^n by a rounded product; the exponent's
point offset dropped; the shift ignored; every point counted as a leading digit: 9 of 9 caught.

*TVB-03, the Bernoulli function (§11.8, TB-Q9; `include/NiTCAD/verify/reference/bernoulli.hpp`,
`tests/verify/benchmark/bernoulli_test.cpp`).* Source read before code: OEIS A000367 and A002445 b-files
(B₂ₖ numerators and denominators, n = k, offset 0 = B₀; A000367's title "Numerators of Bernoulli numbers B_2n");
B₂…B₃₄ typed (every numerator exact in a double). Direct method x / expm1(x) in double-double for |x| ≤ 670;
above, the scaled form (x e^−r)·2^−m, r = x − m ln 2, so B is held below the smallest double; below −670, −x.
Second method: the series 1 − x/2 + Σ B₂ₖ x^2k / (2k)!, k ≤ 17, |x| ≤ 1/2. Values are returned as a double-double
times 2^e and rounded only in the comparison. Contract: E, I3, PA, G3, applicability "B(x) a normal double".
Choice: the TB-Q4 rule (reference uncertainty ≤ 1e-3 of the criterion) applied to TVB-03's self-checks too.

Measured (Release), every outcome as frozen. Reference self-checks: series vs direct 0.88 u_dd, B(−x) = B(x) + x
0.31 u_dd, scaled vs direct on [600, 670] 188 u_dd (r's absolute error ~|x| u_dd from m ln 2 becomes e^−r's relative
error; the direct exp carries ~|x|/8 u_dd too), all ≪ the 8.9e-19 budget; B(1) = 1/(e − 1) to 0. Code: 12 424 points
judged, worst 0.99 ε at x = −0.794 (criterion 4 ε); 602 points refused as subnormal B (diagnostic, worst relative 1:
the code returns 0). **Finding (iii), confirmed:** 104 fails on [709.80, 714.95], where std::expm1 overflows and
`assemble::bernoulli` returns 0 while B is a normal double (outcome fail, as frozen; hidden `[.finding]` case fails
with 104). x/ΔV_T ≈ 710 is an 18 V drop on one edge at 300 K; the fix (the scaled form, or x e^−x for x > 709) is
for the owner (VR9). B′ diagnostic: worst 3.2e-14 relative (146 ε) at x = 1.12e-2, just above the series switch,
consistent with the header's "about 1e-13" (the subtraction 1 − x − B).
Mutations (14, built, run, restored): reference — B₁₂'s numerator, B₆'s denominator, the series' −x/2 sign, the
scaled exponent, the negative tail, the normal floor, the comparison's scaling, the factorial: caught; B₂₄'s
numerator and B₃₄'s denominator survive (their terms are below u_dd at |x| ≤ 1/2, under the budget by ~10¹³), the
B′ series order survives (diagnostic, not asserted). Code — x/(e^x − 1) for expm1, 1/12 → 1/12.1: caught; the
series switch moved to 0.1 survives (the series is still within 4 ε there: equivalent under the criterion).
11 of 14 caught, 3 explained.

*TVB-02, Fermi–Dirac integrals and the inverse (§11.8, TB-Q4–Q8; `include/NiTCAD/verify/reference/fermi_dirac.hpp`,
`src/verify/reference/fermi_dirac.cpp`, `tests/verify/benchmark/fermi_dirac_test.cpp`).* Sources read before code:
DLMF §25.12 (v1.2.8): 25.12.14 (definition, 1/Γ(s+1), s > −1), 25.12.16 (F_s = −Li_{s+1}(−e^x)), 25.12.10 (Li
series); OEIS A078434 (ζ(3/2)), A059750 (digits of −ζ(1/2)), A247041 (ζ(5/2)), 40 digits each (A059750 and A247041
state OFFSET 1,2 for a 1.xxx value; the point placed as their EXAMPLE lines show). Form: t = u², cut at
U = √(η + 100), split at √η for η > 0; tail bounded by integration by parts. **Correction during build:** the first
draft's bound used the p ≥ 1 recurrence for p = 0 (F₋½), where the coefficient (p − 1)/2 is negative and the
inequality reverses; the added test that integrates the omitted tail [U, √(η + 200)] caught it (actual tail above the
bound at η = −40); the bound is now e^−U²/(2U) for p = 0. Primary tanh-sinh, second method adaptive Gauss–Kronrod,
epsrel 1e-27 a piece; U_ref = max(|TS − GK|, TS estimate) + tail bound. Inverse: Brent on F½,ref(η) = ν in
double-double, bracket ±1e-10·max(1, |η|) around the code's answer (widened ×100 until it brackets), t = 1e-18·scale.
Contract: E, I3, PA, G4, second method named. Grid: 320 points (the ruling's "321" counted 40.0123, outside the
range) plus −2 ± 1e-9, 40 − 1e-9, 0.

Measured (Release). Reference: uncertainty ≤ 1e-15 relative at every point (the budget); series for η ≤ −2 and ζ
anchors (1.8, 0.6, 0.1 u_dd) agree; tail bound holds at η = −40, −2, 0, 20, 40 for each order; TVD-14, recorded: the
existing gate's helper (`tests/physics/references.hpp`) agrees to 1.0e-16 (F½) and 9.4e-17 (F₋½) everywhere.
Code: F½ worst 9.8e-15, F₃/₂ 6.4e-15, inverse 8.4e-15 (1296 Brent evaluations): pass. **Finding (ii), widened:**
F₋½ (`fermi_half().derivative`) misses at 82 points from η = 13.76, 2.1e-12 worst up to 15, 2.0e-11 at 38.5; the
header's "1.3e-13 up to 15, 3e-12 up to 40" and the existing gate's 2e-13/5e-12 hold only at that gate's points. The
grid's offset (η_k ≡ 0.0123 mod 0.025) puts every point at 0.49 of the code's table cell, near where the quintic
Hermite derivative's error peaks; the existing points are scattered. Diagnostic: F₃/₂'s derivative vs F½ 8.4e-12.
Mutations (12, built, run, restored): F½'s normalisation, the cut L = 30, the tail bound zeroed, the p = 0 draft
bound, F₃/₂'s power, the η = 0 factor, the logistic branch, the series' alternation and stop, a 1e-14 integrand
factor, epsrel 1e-10: caught; no split at √η survives (accuracy holds without it: equivalent). 11 of 12.
Ctest `verify_benchmark.fermi_dirac`: Release 2.7 s, Debug 18.6 s. Verify suites 28/28 in both configurations;
serial Debug 49.6 s, within §11.9's 60 s.

TVB-01–03 complete (§11.8). Findings for the owner: (i) ħ truncation, (ii) F₋½ widened, (iii) B(x) for
709.78 < x < 715; none changes code in TVR-1. Next in §11: the report and CLI (§11.7).

*Report and CLI (§11.7, RC-Q1–Q7).* **Benchmarks as library functions (RC-Q1):** `verify/benchmark/registry.hpp`
(`benchmarks()`: TVB-01, 02, 03, 32, all fast; `run_tvb01/02/03/32`), sources `src/verify/benchmark/{registry,
tvb01,tvb02,tvb03,tvb32}.cpp`; `tvb32.hpp`'s `tvb32_entries()` keeps the per-check results for the test. The tests now
call these and hold their outcomes to the frozen expectations: `tests/verify/benchmark/{constants,bernoulli,
fermi_dirac}_test.cpp` and `tests/verify/consistency/decks_test.cpp` (same per-check worst ratios as before the move:
pn C1 0.506, C3 0.0085; Schottky C3 0.0097; MOS-C C8 0.0056; entries 92/176/162). TVB-02's
test reads the library's `q_ref`/`u_ref` instead of recomputing the references (the series and TVD-14 checks compare
in double, the half-ulp within the 1e-15 budget). **Known findings (RC-Q2):** `known_findings()` holds (i) TVB-01
"hbar" everywhere, (ii) TVB-02 "F_-1/2" everywhere, (iii) TVB-03 "B(x)" for x > ln DBL_MAX; `apply_known_findings`
turns a fail in the region into `finding` with the ruling cited and adds the "no longer reproduces" note when a region
holds no fail. **io (RC-Q3):** `io::Content::verify_report = 7`, `include/NiTCAD/io/verify_report.hpp` (plain data:
io does not see the verify layer) and `src/io/verify_report.cpp` (the existing container's `write_file`/`read_file`;
keys local to the kind; text lists as one-field records); existing io code and tests unchanged. **Report
(RC-Q4):** `make_report` (identity and grades as "class independence provenance; citation; version", U_num "not
estimated", the two build-time evidence lines) and `report_csv` (one row per comparison, RFC-4180 quoting, numbers by
`std::to_chars`); threads 1 (the deck solves use the default Eigen backend). **Decks (RC-Q5):** CMake writes
`build/<cfg>/generated/verify_decks.hpp` from the three decks (raw strings, `file(CONFIGURE)` so an unchanged deck
does not rebuild; the decks are configure dependencies). **CLI (RC-Q6):** `src/cli/verify.{hpp,cpp}`, dispatched from
`main.cpp` (usage line and header added); a line per benchmark (overall outcome, counts, notes), the summary line,
`--report <file>` writes the report and `<file>.csv`. **VR7/VR8 (RC-Q7):** `tests/cli/verify_probe.cpp` (the shipped
registry plus `TVB-TEST-FAIL`; `compare a b` with the clock cleared) and `tests/cli/verify_exit.cmake`; ctests
`verify_cli.{finding_is_0, unknown_id, unknown_tier, unknown_option, unwritable_report, injected_fail,
reproducible_tvb32}`; in-process `tests/verify/benchmark/registry_test.cpp` (6 cases: registry, known-finding region
and other-benchmark isolation, stale note, overall and summary, CSV quoting, io round trip and VR7 with TVB-01/03).

Measured. `nitcad verify` (Release): TVB-01 finding (12 comparisons; 1 finding, 11 pass), TVB-02 finding (1296; 82
finding, 1214 pass), TVB-03 finding (13 026; 602 NA, 104 finding, 12 320 pass), TVB-32 pass (430; 220 NA, 210 pass);
"verify: 4 benchmarks; 1 pass, 3 finding, 0 fail, ..."; exit 0; CSV 14 764 rows. Shown to detect what they claim:
`verify` exiting 0 on a fail fails `verify_cli.injected_fail`; a report field that varies run to run fails
`verify_cli.reproducible_tvb32`. Verify ctests 36/36 in Release (4.8 s) and Debug. **VR10 miss, for the owner:**
the serial Debug verify run is 63.2 s against §11.9's 60 s; `verify_cli.reproducible_tvb32` (two TVB-32 runs) takes
11.7 s of it (others: `verify_benchmark.fermi_dirac` 18.7 s, `verify_consistency.transient` 13.1 s). The budget is
not moved (TVM2).

*VR10 after the ruling (option A).* `verify_cli.reproducible_tvb32` replaced by `verify_cli.reproducible_TVB01` and
`_TVB03` (Debug 0.06 s and 2.2 s). Verify ctests 37/37 in Release (4.9 s) and Debug; serial Debug verify run 57.4 s,
within 60 s. The serial Debug time varies by a few seconds from run to run (63.2 s and 66.7 s were measured for the
previous set), so the margin is about 3 s. Full suites (parallel, -j 6): Release 166/166 (823 s), Debug 155/155
(1091 s). VR1–VR10 met; TVR-1 complete. Open for the owner: findings (i)–(iii) (no code changed in TVR-1, VR9).


#### 12. Findings (i)–(iii): causes and fix proposals (proposed 2026-11-05; rulings before code)

TVR-1 changed no code under test (VR9). Each fix below is its own code change, so these rulings come before any of
it. Proposed vehicle (FQ0): one small unit, Unit 56, on `physics/reference-findings`, branched from
`verify/reference-numerics` so that TVB-01/02/03 judge the fixes. Each fixed finding's known-findings entry is
removed in the same commit. RC-Q2 makes a stale entry report "no longer reproduces", so a missed removal shows up.
The frozen criteria do not move (TVM2).

**Finding (i): ħ off by 6.13e-10 (TVB-01).**
- Cause: `include/NiTCAD/base/constants.hpp:21` stores `1.054571817e-34`. That is CODATA 2018's printed ħ, rounded to
  10 significant figures. ħ = h/2π with h = 6.62607015e-34 exact, so ħ = 1.054 571 817 646 156 5…e-34. The
  comment "(exact)" is wrong: the value is exact, the stored decimal is not.
- FQ1 (recommended A):
  - A: replace the literal with the correctly rounded double `1.0545718176461565e-34`, with a comment giving
    h/2π. TVB-01 checks it to its γ_n bound.
  - B: compute `h_J_s / (2π)` in constexpr, adding `h_J_s = 6.62607015e-34`. This can be up to about 1.5 ulp
    off, where A is correctly rounded.
- Effect: the 10 uses (band-to-band, gate tunnelling, MLDA, the thermal γ_ii) move by 6.1e-10 relative. Any
  existing expectation that moves is reported to the owner with its old and new values. None is re-baselined
  silently.

**Finding (ii): F₋½ misses 1e-12 from η ≈ 13.76, worst 2.0e-11 at η = 38.5 (TVB-02).**
- Cause: F₋½ is returned as F½·(1 + d lnγ/dη), with the slope taken from the quintic Hermite interpolant of lnγ.
  - The interpolant's derivative error is O(h⁵), with step h = 0.025.
  - As η grows, 1 + d lnγ/dη ≈ 3/(2η). The slope's absolute error is therefore amplified about η/1.5 times
    (≈ 26 at η = 38.5) into F₋½'s relative error.
  - TVB-02's grid sits at 0.49 of a table cell, where the interpolation error is largest. The existing gate
    (`tests/physics/fermi_dirac_test.cpp:79`, 5e-12 above η = 15) samples elsewhere in the cell. So the
    header's "3e-12 up to 40" is optimistic.
- FQ2 (recommended A):
  - A: reduce `step` in `src/physics/fermi_dirac.cpp` from 0.025 to 0.01. This applies to both the F½ and F₃/₂
    tables, which share the constant.
    - Expected slope error is about 2.5⁵ ≈ 98 times smaller, so the worst case becomes about 2e-13: a 5× margin
      on 1e-12.
    - Each table grows from 1681 to 4201 nodes (about 100 KB each).
    - The one-time build cost per process grows about 2.5×. It will be measured and reported. If it is
      material for a short run, that comes back to the owner.
  - B: step 0.0125, about 32× smaller. The margin is only about 1.6×, too thin.
  - C: interpolate F₋½ separately. Rejected: it breaks the Unit 14 property that the returned derivative is
    exactly that of the returned value, which is what the Newton Jacobian relies on.
- FQ3 (recommended): correct the accuracy statements in `include/NiTCAD/physics/fermi_dirac.hpp` to the values
  measured after the change. Leave `fermi_dirac_test.cpp`'s bounds as they are (they only get looser than the new
  accuracy). TVB-02 is the mid-cell gate.
- Effect: F½ and F₃/₂ values in (-2, 40] change at rounding level (~1e-15). Fermi-Dirac results can move in their
  last digits. Moved expectations are handled as in FQ1.

**Finding (iii): B(x) = 0 on x ∈ (709.78, 714.95] (TVB-03, 104 fails).**
- Cause: `include/NiTCAD/assemble/bernoulli.hpp:30` returns `x / std::expm1(x)`. `expm1` overflows to inf above
  ln DBL_MAX = 709.78, so B returns 0 while x·e⁻ˣ is still a normal double. The header's "B(x) = 0 once x e^-x
  underflows" is wrong on that interval. `bernoulli_derivative` has the same defect: there it should return
  (1 − x)e⁻ˣ. TVB-03 does not compare B′.
- FQ4 (recommended A):
  - A: for x > 709, return B = (x·e^(−x/2))·e^(−x/2) and B′ = ((1 − x)·e^(−x/2))·e^(−x/2).
    - x/2 is exact, and each factor stays normal.
    - The bound is about 2 ulp from `exp` plus 1 ulp from the products. That is under the 4 eps criterion; the
      measured worst case will be reported.
    - Below 709 the code path is unchanged, so every existing result stays bit-identical.
  - B: use the same form for all x > 40. It is equally accurate, but it changes values at ulp level in
    (40, 709] for no gain.
  - C: accept the finding as a limitation. In the Scharfetter-Gummel flux, x > 709 means a 709 kT (18 V at 300 K)
    drop across one edge, where both 0 and 1e-308 are effectively zero.
- FQ5 (recommended): add a B′ case to the unit's tests on (709, 745], against a double-double reference of B′ built
  like TVB-03's scaled form. TVB-03 itself stays frozen.

**Order and acceptance (FQ6, recommended):**
- Order: (iii), then (i), then (ii), one commit each.
- Each commit:
  - removes its known-findings entry;
  - shows the benchmark passing;
  - runs the ctests of the touched areas.
- Both full suites run before the unit is offered.
- Nothing is pushed without the owner's word.
#### 13. Unit 56: the fixes of findings (i)–(iii) (on `physics/reference-findings`)

**Rulings:** FQ0–FQ6 are approved as recommended (owner, "go", 2026-11-05, after "what do you recommend"). In the
same instruction the owner authorised pushing `verify/reference-numerics` and opening the TVR-1 PR (#76). Unit 56 is
stacked on that branch and will be rebased onto main once #76 merges.

*As built, finding (iii) (FQ4, FQ5):*
- `assemble::bernoulli` and `bernoulli_derivative` above x = 709 return `(x·h)·h` and `((1 − x)·h)·h`, with
  h = e^(−x/2). The code path below 709 is unchanged, so every value there stays bit-identical.
- TVB-03 now passes at all 12424 judged points, worst 1.19 eps at x = 713.6. The 602 subnormal points remain
  not applicable.
- The new B′ case on (ln DBL_MAX, 745] judges 104 points, worst 1.20 eps at x = 713.8.
- Known-findings entry (iii) is removed, and the table has 2 entries.
- The frozen TVB-03 expectation in `tests/verify/benchmark/bernoulli_test.cpp` changes from "finding where expm1
  overflows" to "pass everywhere". The hidden `[.finding]` case is gone: it is now the visible case.
- `registry_test` expects TVB-03 to report `pass` and the table to hold 2 entries.
- Ctests: verify 37/37, assemble 18/18 (Release).

*As built, finding (i) (FQ1):*
- `base::hbar_J_s` is now `1.0545718176461565e-34`, the correctly rounded h/2π.
  - TVB-01 passes 12/12. The new case checks that it is the double nearest the double-double h/2π.
  - Known-findings entry (i) is removed, leaving 1 entry.
- Existing expectations that moved; each is listed here with old and new values:
  - `tests/base/constants_test.cpp:13` pinned the old literal; it now pins the new one.
  - `tests/physics/quantum_confinement_test.cpp` QC2 wrote out the old ħ in its reference; it now uses h/2π.
  - Three cases sit at the rounding floor. They were ruled by the owner on 2026-11-05, choosing the recommended
    option each time:
    - `tests/solve/quantum_confinement_test.cpp` QC5 (reported, not gated, but a solve that does not converge
      fails it): with the model on at V_G = −2 V, the equilibrium solve no longer reaches tol_update 1e-13
      within 100 iterations, though it does with the old ħ. It converges at 2e-13, so tol_update is raised from
      1e-13 to 2e-13.
    - `tests/deck/gate_tunnelling_test.cpp` GC7: the last correction is 3.77e-13, against 3.24e-14 with the old
      ħ. The window cutoff is raised from 1e-13 to the sweep's tol_update of 1e-12.
    - GC11: at 0.25 V the 1D/2D agreement is 1.24e-9, against 7.84e-10 with the old ħ. The bound goes from 1e-9
      to 2e-9. tol_update 1e-12 does not converge in 2D.
  - A 6e-10 change in ħ alters GC7's Newton path visibly (iteration 2: 1.03e-2 → 1.61e-2). The gate-tunnelling
    capacitor is near-singular (GR5), which is the likely reason. Recorded, not investigated.
- Release suite with (i) and (iii): 163/166. The three failures were base, solve_mlda and gate_tunnelling,
  before the edits above; after them all three pass, as do physics 2/2 and verify 37/37.

*Finding (ii): FQ2 falsified; accepted as a limitation (owner, 2026-11-05, the recommended option, after measurement).*
- Step 0.01 made F₋½ worse: worst 2.28e-11 with 91 misses, against 2.0e-11 with 82 at step 0.025.
- The error is not interpolation. F₋½ is the slope of the tabulated ln F½, and the slope term d0·(y0 − y1)/step
  turns the quadratures' relative error (about 3e-15 per node) into F₋½'s relative error, scaled by about
  η/(1.5·step).
- A prototype stored per-cell differences ln(f₍ᵢ₊₁₎/fᵢ) − Δη instead of the rounded ln γ, to avoid ln γ's
  rounding at its magnitude (|ln γ| ≈ η). It halved the worst error to 1.1e-11 but raised the misses to 104.
  Reverted, not committed.
- Meeting 1e-12 would need a table built above double precision. That is a new design and was not taken.
- Where F₋½ is used: only as a derivative, in the Jacobian entries of `physics/statistics.cpp` and in the
  energy-balance degenerate terms (`assemble/energy_balance.cpp:108`). It exceeds 1e-12 only above η ≈ 14
  (n ≳ 1e21 cm⁻³ in silicon).
- Done:
  - FQ3: the accuracy statement in `physics/fermi_dirac.hpp` now gives the mid-cell values measured by TVB-02:
    7e-14 up to η = 5, 2.1e-12 up to 15, 4.3e-12 up to 25, 2.0e-11 up to 40. F½ is within 1e-14.
  - The text of known-finding entry (ii) records the acceptance.
  - Physics code is unchanged. TVB-02 and its tests are unchanged, and the hidden `[.finding]` case still fails
    by design.

*Unit 56 close-out.* Full suites on e5e9847, run in parallel (-j 6): Release 166/166 (805 s), Debug 155/155 (1070 s).
In `nitcad verify`, TVB-01 and TVB-03 pass, TVB-02 reports the accepted finding (ii), and TVB-32 is unchanged.
The unit is stacked on PR #76. It will be rebased onto main once #76 merges, then offered as its own PR on the
owner's word.


#### 14. TVR-2 (Unit 55b) proposal: code verification core (proposed 2026-11-05; on `verify/code-verification`; rulings before code)

TVR-2 follows the roadmap in §9 and catalogue §6: TVB-04, 05, 06, 07, 09, 11, 12 and 14, the TVQ3 hook, and the
`verify/study` driver. This section proposes how to build them. It also reports one place where the plan does not
fit the code (CV-Q3). The criteria (§4.4) are frozen in a criteria commit after these rulings and before any code
(TVM2).

**14.1 Scope.** In:
- (a) `verify/study`: mesh ladders per family, the Newton-tolerance check, and GCI driving. The formulas are
  TVR-1's `verify/reference/metrics`.
- (b) The TVQ3 manufactured-source hook in `assemble` (approved with TVQ1–TVQ6 on 2026-10-31).
- (c) The references, each with L0 tests:
  - manufactured solutions and their sources;
  - the equilibrium pn junction first integral (TVB-07);
  - the Kingston–Neustadter MOS charge (TVB-09);
  - the resistor (TVB-11);
  - the constant-field exponential (TVB-12).
- (d) The benchmarks, registered in TVR-1's registry with tiers, and their tests.

Out:
- TVB-05 and 06, if CV-Q3 rules B.
- Unit 17's use of the studies (TVR-6).
- Any change to the existing tests, gates, tolerances or solver beyond (b) and a ruled CV-Q3.

**14.2 A finding before code: TVQ3's premise is false.** TVQ3 states that the hook "is the only change to existing
code the framework needs for its P1 set". That is not so:
- `src/device/device.cpp:75` refuses a device with no semiconductor region ("device needs at least one
  semiconductor region").
- `Device::Device` (`device.cpp:472`) chooses its reference node assuming a semiconductor node exists.
- So the capacitors of TVB-05 (coaxial, spherical) and TVB-06 (dielectric stack), which are pure dielectrics
  between electrodes, cannot be built.
- TVB-04's "linear Poisson" is affected for the same reason.

CV-Q2 and CV-Q3 rule on this.

**14.3 The hook (TVQ3; design for ruling CV-Q4).**
- `PhysicsModels` gains `std::shared_ptr<const ManufacturedSources> manufactured`, null by default.
  - It holds per-node physical sources: a charge density for the Poisson row [C/cm³], and electron and hole
    generation for the continuity rows [cm⁻³ s⁻¹].
  - The sources do not depend on the state, so they add no Jacobian entries.
- It enters as one node term, `manufactured source`, after `storage` in the Unit 40 composition order, and is
  active only when the pointer is set.
- When absent, the composition is unchanged, so every existing result stays bit-identical. A test asserts that F
  and J are bit-equal with the pointer null and with no term.
- It is not in the deck vocabulary, and the header documents it as verification-only.
- The equilibrium-Poisson assembler gets no hook: every MMS case runs the drift-diffusion system (14.4).

**14.4 The benchmarks.** Each has its reference, mesh families, norms and p_th per §5.3. Tolerances are frozen in the
criteria commit.

| ID | Construction | Reference | Families (fast / full tier) | Judged |
|---|---|---|---|---|
| TVB-04 | Poisson MMS. Per CV-Q2 A: Poisson–Boltzmann with a silicon/oxide ε jump on a planar interface. ψ_exact = ψ_eq(N) + a smooth bump vanishing at the contacts. The Poisson source is the exact charge residual, and the continuity rows hold to rounding at equilibrium. | Analytic, with derivatives checked by complex step | 1D tensor / 2D and 3D tensor, 2D generated, 3D generated, swept | L∞ and L2 (p_th 2, tensor); H1 (p_th 1, unstructured); L2 on unstructured families reported only |
| TVB-05, 06 | Capacitors (CV-Q3) | C = 2πεL/ln(b/a); 4πε ab/(b−a); series ε | — | per CV-Q3 |
| TVB-07 | Equilibrium abrupt pn junction, Boltzmann, 1D, with the junction on a dual-cell face | First integral: E(ψ) in closed form, x(ψ) by tanh–sinh, ψ_j by Brent; Chebyshev collocation as the I3 cross-check | 1D tensor ladder (fast) | V_bi to rounding (D); E_max and ψ(x) in L∞, p_th 2; depletion edges reported |
| TVB-09 | MOS-C, uniform substrate, 1D, Boltzmann, lumped gate | Kingston–Neustadter F(ψ_s) in closed form; V_G = V_FB + ψ_s − Q_s/C_ox; C_LF by complex step | 1D tensor ladder | Q_s(V_G) and C_LF relative error, with p̂; V_FB to rounding |
| TVB-11 | Uniform resistor, ohmic contacts | I = q μ N A V / L | 1D and 2D tensor | D: to c·ε, c = 64 |
| TVB-12 | SG constant-field exactness: 1D, no G−R; the hook's Poisson source makes ψ linear and n exponential | The exact discrete solution | 1D, several Péclet numbers | D: nodal n to c·ε |
| TVB-14 | Full drift-diffusion MMS: ψ, φ_n and φ_p smooth; Boltzmann, SRH on; the three sources taken from the exact state | Analytic; sources in double-double, derivatives checked by complex step | 1D and 2D tensor / 2D and 3D tensor, 2D generated, 3D generated, swept | ψ, φ_n and φ_p in L2 and L∞ (p_th 2, tensor) and H1 (p_th 1, unstructured) |

**14.5 Owner questions (recommendations in the last column).**

| # | Question | Recommendation |
|---|---|---|
| CV-Q1 | Order of work: hook, then study driver, then TVB-12 and 11, then 04, then 07 and 09, then 14. One commit each, with its references' L0 tests first? | Yes |
| CV-Q2 | TVB-04's form: (A) a Poisson–Boltzmann MMS with an oxide region (the semiconductor path; no device change), or (B) the linear Poisson the catalogue names, which needs CV-Q3 A? | A. It exercises the same operator and ε jump, with the carriers' charge exact in the source |
| CV-Q3 | Capacitors TVB-05 and 06. (A) A narrow device change: allow insulator-only devices whose every component has an electrode, with the reference node chosen among electrode nodes. This touches existing code; it gets its own tests and leaves current devices bit-identical. (B) Defer both to TVR-4 and record the dependency, as TVD-15 does. (C) Model the plates as degenerate semiconductor: rejected, not exact | B. TVB-04 already covers the Poisson operator and the unstructured families, and B keeps TVQ3 the only change to existing code |
| CV-Q4 | The hook as designed in 14.3: per-node, state-independent sources, one node term after storage, null by default, not in the deck | As 14.3 |
| CV-Q5 | Source discretisation: the scheme's own lumping S(x_i)·V_i (the way doping charge enters), or dual-cell quadrature? | Lumping. It is the scheme's own consistency, O(h²) on a smooth S |
| CV-Q6 | Orders: use §5.3's working assumptions (p_th 2 for L2/L∞ on tensor, 1 for H1 on unstructured) and gate them per §4.4 E, while the Eymard–Gallouët–Herbin theorem is still unread? | Yes. The report marks the citation [unverified]. A miss is a fail, and goes to the owner before any criterion moves |
| CV-Q7 | Ladders, r = 2: 4 levels for tensor in 1D, 2D and 3D; 4 levels for 2D generated; 3 levels for 3D generated and swept (full tier). Unstructured h = (V/N)^(1/D) | Yes |
| CV-Q8 | Iterative error: U_iter is the final Newton correction propagated to the metric. It must be ≤ U_disc/10 before an order is computed (§5.2). A rounding-floor stall (as in Unit 56) goes to the owner, not to a looser tolerance | Yes |
| CV-Q9 | Tiers and time: TVR-1's fast tier already takes 57.4 s of the 60 s serial Debug budget (VR10). TVR-2 goes in `full`, plus one coarse case per benchmark in `fast`. The fast tier's budget is kept by moving nothing else in. The full tier gets its own ctest label and the existing ctest timeouts | Yes |
| CV-Q10 | Criteria freeze: the tolerance on the finest mesh per benchmark and the D-class constant c, in a separate commit after these rulings and before code | Yes |
| CV-Q11 | Unit numbering and branch: Unit 55b, on `verify/code-verification` | Yes |

**14.6 Rulings (owner, "choose recommended", 2026-11-05).** CV-Q1 through CV-Q11 are approved as recommended:
- TVB-04 is a Poisson–Boltzmann MMS with an oxide region.
- TVB-05 and 06 are deferred to TVR-4. Their dependency, insulator-only devices (`device.cpp:75`, `:472`), is
  recorded as TVD-15 does its dependencies.
- The hook is built as 14.3 describes, and sources are lumped.
- The §5.3 orders are gated, with the theorem's citation marked [unverified].
- The ladders are 14.5's CV-Q7, and the iterative-error rule is CV-Q8.
- The fast/full split is CV-Q9.
- The criteria are frozen below (14.7) before any code.

**14.7 Criteria (frozen 2026-11-05, before code; TVM2: never widened after a run).** Notation:
- ε = 2^-52.
- N is the number of nodes along the current path (1D) or per axis (tensor).
- Relative norms are normalised by the same norm of the exact solution. For ψ the normaliser is the potential
  scale: V_bi for TVB-07, and max|ψ_exact − ψ_eq| + V_T for TVB-04 and TVB-14.
- p̂ comes from §4.2's observed order on successive pairs. The asymptotic-range check of TVR-1 applies with ≥ 4
  levels; with 3 levels it is `not_evaluated` and p̂ is judged on the finest pair.
- No order is computed unless U_iter ≤ U_disc/10 (CV-Q8). If that cannot be met, the outcome is `not_evaluated`
  and the case is brought to the owner.

| ID | D criterion | E (a) order | E (b) accuracy on the finest level | Ladder (fast / full) |
|---|---|---|---|---|
| TVB-11 | relative \|I − I_ref\| and nodal ψ in L∞, each ≤ 64 ε N² (N²: the 1D Laplacian's condition) | — | — | 1D N = 32, 128 / + 2D 32² |
| TVB-12 | nodal n relative ≤ 64 ε N², at Péclet numbers \|ΔΨ/V_T\| per edge ∈ {0.01, 1, 10, 40} | — | — | 1D N = 64 (fast) |
| TVB-04 | — | tensor L2 and L∞: p̂ ∈ [1.8, 2.3]; unstructured H1: p̂ ∈ [0.8, 1.3]; unstructured L2 reported | tensor L2 relative ≤ 1e-4; unstructured H1 relative ≤ 1e-2 | 1D N = 32, 64, 128, 256 / + 2D tensor 16²…128², 3D tensor 8³…64³, 2D generated 4 levels, 3D generated and swept 3 levels |
| TVB-07 | V_bi (contact ψ difference) relative ≤ 64 ε | ψ L∞ and E_max: p̂ ∈ [1.8, 2.3] | ψ L∞ relative ≤ 1e-4; E_max relative ≤ 1e-3 | 1D N = 64, 128, 256, 512 (fast: the coarsest two, order not judged) |
| TVB-09 | V_FB relative ≤ 64 ε | Q_s(V_G): p̂ ∈ [1.8, 2.3] at each V_G | Q_s relative ≤ 1e-3; C_LF relative ≤ 1e-3; on a V_G grid from accumulation to strong inversion, with \|Q_s\| ≥ 1e-3 of its range (below that, `not_applicable`, predicate named) | 1D N = 64…512 (fast: one bias) |
| TVB-14 | — | as TVB-04, for ψ, φ_n and φ_p | as TVB-04, for each of the three | as TVB-04 |

References are analytic (U_ref at rounding) or computed in double-double, with U_ref ≤ tolerance/10 checked by
§4.4 E (c). The TB-Q4 rule applies: a reference whose uncertainty exceeds 1e-3 of the criterion is
`not_evaluated`.

**14.8 As built.**

*Step 1, the hook (CV-Q4).*
- `assemble::ManufacturedSources` and `PhysicsModels::manufactured` (`models.hpp`); the node term `manufactured source` after `storage` (`terms.cpp`).
- `DriftDiffusion::create` refuses a source vector of the wrong size or with non-finite values, and refuses the hook with electrothermal or energy balance.
- `tests/assemble/manufactured_test.cpp` (4 cases) checks:
  - absent, the composition, F and J are equal to the plain system's;
  - present, the term alone writes +V S/Ns, +V G_n/R0 and −V G_p/R0 within 1e-14, skips contact and electrode rows, puts charge on interior oxide nodes, and adds no Jacobian values;
  - an empty vector leaves its rows bit-equal;
  - the refusals.
- Measuring the term as a difference of two full residuals is wrong: the later flux terms cancel, so their rounding shows (2.2e-14 relative was seen). The test therefore evaluates the term alone.
- assemble tests: 154/154 pass (Release).

*Step 2, `verify/study` (CV-Q6–Q8).*
- What it is: `nitcad_verify_study`, which depends only on `mesh` and `verify_reference`; the benchmarks run the solves.
- What it computes:
  - `order_study`: the observed orders from errors against the exact solution, with p̂ from the finest pair.
  - The CV-Q8 rule: U_iter ≤ e/10 at every level, else the order is not evaluated.
  - The asymptotic check from 4 levels.
  - `judge_order` against [p_th − 0.2, p_th + 0.3].
  - `mean_spacing` = (V/N)^(1/D).
- **Interpretation for the owner (14.7, "TVR-1's asymptotic check applies").** TVR-1's check has two tests: the order spread (≤ 0.1), and the GCI ratio ρ_k ≈ 1.
  - The GCI ratio is relative to a nonzero solution value. Applied to errors that converge to 0, it gives ρ_k = r^−p (0.25 at second order), so every ladder would fail it.
  - For studies against an exact solution, the check is therefore the order spread alone.
  - This is narrower than TVR-1's two-test form. It is recorded here and not hidden.
- Tests: `tests/verify/study/study_test.cpp`, 7 cases. They cover exact orders at r = 2 and r = 1.7, the band edges, 3 levels without the check, a drifting ladder that is not judged, the CV-Q8 boundary, invalid input, and the mean spacing in 1D, 2D and 3D.

*Step 3, TVB-12 and TVB-11 (D class).*
- **Change to step 1 (14.3 said `shared_ptr`).** `PhysicsModels::manufactured` is a non-owning `const ManufacturedSources*`. A `shared_ptr` member made `PhysicsModels` a non-literal type and broke the existing `constexpr PhysicsModels` uses (`tests/solve/heterojunction_test.cpp:87`, `tests/assemble/energy_balance_test.cpp:497`). The caller keeps the sources alive through the solve; the hook tests and the benchmarks own them and pass `.get()`.
- `src/verify/benchmark/tvb12.cpp`, as frozen in 14.4/14.7:
  - 1D, 64 nodes on [0, 1 µm]; N_D 1e18 at the left contact node, 1e15 at the right, 1e16 inside; constant mobilities, no G−R, no BGN, Boltzmann.
  - The contact values come from an equilibrium solve. The bias puts |Δψ/V_T| per edge at 0.01, 1, 10 and 40. The exact n and p are written without overflow, and the TVQ3 charge n − p − (N_D − N_A) makes ψ linear the discrete Poisson solution.
  - Each solve starts from the exact state perturbed by 1e-3.
  - Result (Release): every nodal n and p error is between 9e-16 and 3.4e-14, against 64εN² = 5.8e-11. 8/8 pass.
- `src/verify/benchmark/tvb11.cpp`: uniform N_D = 1e16, V = 0.1 V, 1D N = 32 and 128 and 2D 32 × 32 (per unit width, width 1 µm); μ the material's mu_max. Besides the 14.7 criterion, the current comparison also requires the two contact currents to have opposite signs.
  - Result (Release): current errors 1.2e-14 to 2.3e-14, ψ-linearity about 1e-16; criteria 1.46e-11 (N = 32) to 2.33e-10 (N = 128). 6/6 pass.
- Registry: 6 shipped benchmarks; TVB-11 and TVB-12 are `fast`. CV-Q9 put TVR-2 in `full` with one coarse case per benchmark in `fast`. These two are already single, coarse D-class runs (26 s is the whole serial Debug `nitcad verify`), so each is its own fast case. `registry_test` now expects 6 ids (was 4); this is the registry's own list, not a moved expectation of a benchmark.
- `tests/verify/benchmark/code_verification_test.cpp` (`[tvr2]`): both benchmarks, every comparison pass (frozen outcomes).
- Time: the verify ctests take 39/39 in Release. Serial Debug takes 58.2 s against the 60 s budget (VR10; it was 57.4 s), so the fast tier has under 2 s left. TVB-04, 07, 09 and 14 go in `full` only, per CV-Q9.

*Step 4, TVB-04 (tensor families) and tier selection.*
- **Tiers.** `Benchmark::run` takes the requested tier, and `Benchmark::tier` is the lowest tier a benchmark runs in. `nitcad verify --tier t` runs every benchmark whose tier is at or below t. With `--id`, a benchmark above t is refused with exit 64. The TVR-1 and TVB-11/12 entries ignore the tier.
  - CLI tests: `verify_cli.below_tier` and `verify_cli.at_tier`, on the probe's test-only `TVB-TEST-FULL`.
- **Construction** (`src/verify/benchmark/tvb04.cpp`), as 14.4 with CV-Q2 A:
  - SiO2 on [−L, 0] with a metal electrode at −L, and p-Si (N_A = 1e17) on [0, L] with an ohmic contact at L; L = 50 nm.
  - The interface lies at the midpoint of the edges that join the regions (`straddle_interface`, as the decks do it).
  - ψ = ψ_c + A g(x/L) h(y/L) h(z/L), with A = 0.2 V. g = (1−u)e^{2u} in Si and (1+u)e^{βu} in SiO2, with β = (ε_si/ε_ox)(α−1) − 1, so D·n is continuous. h = 1 + cos(πv)/2, which has zero normal flux on the sides.
  - Equilibrium Boltzmann carriers make every SG flux vanish. The TVQ3 charge is S = −εΔψ/q + n − p − N.
  - The electrode is biased so that its potential equals ψ_c.
  - The analytic g′, g″, h′ and h″ are checked against complex-step derivatives on every run (the second method).
- **Tier placement.** fast: 1D N = 32…256. full (and nightly): + 2D tensor 16²…128² and 3D tensor 8³…64³. Per the 14.7 ladder column and CV-Q13 (14.9), this is in full despite TVD-13's "nightly: 3D ladders".
- **3D solves (CV-Q15).** The 3D levels use MKL PARDISO. With Eigen SparseLU, 64³ (786k unknowns) held 26 GB and had not finished after 15 minutes, so that run was stopped.
  - `nitcad verify --mkl <mkl_rt.dll> --threads <n>` loads the runtime and sets PARDISO's threads; the report's `threads` field records n when MKL is loaded.
  - Without MKL, the 3D family's three comparisons are `not_evaluated`, with the missing input named.
  - Side effect: without MKL, TVR-1's stale-finding check prints "known finding (iv) no longer reproduces" for the 3D entry, because its region holds no fail. The note is wrong in that case and is left as is.
- **Full-tier test (CV-Q16).** `tests/verify/benchmark/code_verification_full_test.cpp` (`verify_full.tvb04`, label `full`, run serially) exists only with `NITCAD_VERIFY_FULL_TIER=ON` in Release, and takes the runtime from `NITCAD_MKL_RUNTIME`. CI does not set the option. Run it with `ctest -L full`.
- **Results (Release).**
  - 1D: L2 p̂ = 2.001 and L∞ p̂ = 1.989; order spreads 5e-4 and 0.024; finest L2 7.1e-5 ≤ 1e-4. Pass.
  - 2D tensor: L2 p̂ = 2.002 and L∞ p̂ = 1.972; spreads 0.008 and 0.071. Both orders pass. **Finest L2 (128²) is 2.88e-4 > 1e-4: fail**, brought to the owner (CV-Q12), and now known finding (iv).
  - U_iter/e ≤ 4e-11 at every level.
- 3D tensor (PARDISO, 6 threads):
  - L2: 0.118, 0.0228, 0.00501, 0.00117; orders 1.93, 1.99, 2.00, with p̂ = 2.001 and spread 0.068. Pass.
  - L∞: orders 1.60, 1.84, 1.93; spread 0.23 > 0.1. **Not evaluated** (outside the asymptotic range, CV-Q14). The reason is stated as such, not as a missing input.
  - Finest L2 1.17e-3: known finding (iv).
- `verify_full.tvb04`: every frozen outcome held; overall `not_evaluated`; 710.8 s.
- Time: TVB-04 fast takes 0.17 s in Debug. In Release, 1D + 2D takes 2.1 s, and 64³ takes 709 s at about 20 GB.
- The verify ctests pass in Release and Debug (41/41). Serial Debug takes 59.5 s against the 60 s budget (VR10), so the fast tier is at its limit. The fast cases 14.7 still lists (TVB-07's coarsest two levels, TVB-09's one bias) go to the owner if they take it over 60 s.

*Step 5, TVB-04 (unstructured families).*
- **Meshes** (`src/verify/benchmark/tvb04.cpp`). Each family uses the same domain as the tensor families, [−L, L] × [0, L]^(D−1), with the same manufactured ψ.
  - Patches `x_min` (electrode) and `x_max` (ohmic).
  - Size bound h_max = L/k. Interface resolution s = h_max, the straddle width, as the tensor axis's L/n_side.
  - 2D generated: `generate_mesh(Geometry2D)`, two rectangles that meet at x = 0, k = 16, 32, 64, 128.
  - 3D generated: `generate_mesh(Geometry3D)`, a box with the planar interface x = 0 and one seed per slab, k = 8, 16, 32.
  - 3D swept: `generate_swept`. The (x, y) cross-section (the interface lies in it) is swept along z on uniform lines h_max apart, k = 8, 16, 32.
  - The device takes each generated mesh's P4 certificate. Node regions are the generator's, and are checked against the side of x = 0. A node on x = 0 is refused, since the interface is straddled.
  - h for the order study is (V/N)^(1/D) with V = 2L^D (CV-Q7). The 3D levels use PARDISO (CV-Q15).
- **Norm.** The discrete H1 seminorm is Σ_edges coupling_area (e_i − e_j)²/length (as in `complex_convergence_test.cpp`). It is relative to the same seminorm of the exact ψ. **Interpretation for the owner:** 14.7's potential-scale normaliser (A + V_T) is a voltage, so it is not used for a gradient norm. The seminorm of the exact ψ does not see ψ_c, which is the problem that rule avoids.
- **Comparisons.** Each family records "psi H1 order, …" (p_th 1, band [0.8, 1.3]) and "psi H1 finest, …" (≤ 1e-2). The L2 and L∞ orders are in the notes, marked "reported, not gated". TVB-04's full tier now has 15 comparisons (was 9).
- **Results (Release, PARDISO 6 threads for 3D).**

  | Family | Nodes (coarse → fine) | H1 orders (coarse → fine) | p̂ | Finest H1 | L2 orders (reported) | Outcome |
  |---|---|---|---|---|---|---|
  | 2D generated | 1 087 → 68 406 | 1.21, 1.11, 1.08 (spread 0.099) | 1.08 | 7.9e-4 | 2.29, 1.86, 2.18 | pass, pass |
  | 3D generated | 3 985 → 224 485 | 1.52, 1.29 | 1.29 | 3.3e-3 | 2.05, 2.03 | pass, pass |
  | 3D swept | 2 538 → 143 550 | 1.45, 1.23 | 1.23 | 3.3e-3 | 1.90, 2.25 | pass, pass |

  - U_iter/e ≤ 6.3e-14 at every level.
  - 3 levels, so the asymptotic-range check is `not_evaluated` and p̂ is the finest pair (14.7).
- **Observations, not findings:**
  - The H1 orders fall toward 1 as the mesh is refined, from above. The coarse pairs (1.52, 1.45) lie outside the band. With 3 levels only the finest pair is judged.
  - 3D generated's p̂ = 1.29 is 0.01 below the band's top.
  - 2D's spread of 0.099 is 0.001 under the 0.1 limit.
  - A fourth 3D level would settle both. It is not in the frozen ladder, and it is not added.
- **Full-tier test.** `verify_full.tvb04` expects the 15 outcomes: the tensor ones as before, and the six unstructured ones pass.
  - Time: the 3D generated and swept ladders take about 13 minutes together; `verify_full.tvb04` passes (19 assertions) in 1505 s. It is still local only (CV-Q16).
- The verify ctests pass in Release and Debug (41/41; the fast tier is unchanged). Serial Debug takes 59.6 s against the 60 s budget (VR10).


*Step 6, TVB-07 (equilibrium abrupt pn junction).*
- **Reference** (`include/NiTCAD/verify/reference/pn_junction.hpp`, `src/verify/reference/pn_junction.cpp`, in the reference layer: standard library and base only, VR2).
  - u_n = asinh(N_D/2n_i) and u_p = −asinh(N_A/2n_i), in double-double. **Note:** the catalogue's V_T ln(N_A N_D/n_i²) is their majority-carrier limit. It differs by about (n_i/N)² relative, 7e-14 here, above the 64 ε = 1.4e-14 of the D criterion. So the exact neutral levels are the reference. This is the same V_bi the catalogue names, written exactly, and the criterion is unchanged.
  - First integral on each side: ½u′² = (q/εV_T)[n_b G(t) + p_b G(−t)], with G(t) = e^t − 1 − t (series for |t| < ¼). u_j is found by Brent on Φ_n = Φ_p. E_max = V_T u′(0).
  - x(ψ) uses tanh-sinh in s = ln|u − u_b|, whose integrand tends to the Debye length, so there is no endpoint singularity. Its inverse is a safeguarded Newton in s.
  - **One defect found and fixed before the L0 tests passed:** recomputing the tail as (u_b − t) − u_b cancelled to about 1e-18 relative, so x(s) was noisy and Newton's bracket closed. Φ is now taken from t directly, and a closed bracket is an error, not a silent stop.
  - **I3, a second method:** Chebyshev collocation of the finite-domain problem. It uses 4 + 4 pieces, patching u and u′, damped Newton, and a dense LU in double. It shares no code with the first integral.
  - L0 tests (`tests/verify/reference/pn_junction_test.cpp`, 6 cases, 44 assertions):
    - symmetric doping gives u_j = 0 and an odd profile;
    - the junction condition holds to 6e-32, re-evaluated with G = expm1 − t;
    - bulk neutrality and mass action hold to 1e-30;
    - the profile agrees with the ODE, shot from (u_j, E_max) by dopri5 in double-double on each side to the depletion edge plus 3 L_D, to 1e-15 of V_bi (measured 0 at double precision);
    - E_max agrees with the depletion approximation (with the −2V_T correction) to 1e-3, on junctions depleted on both sides. A 1:100 junction is not, by ~46 %, and is not a case;
    - collocation at degrees 40 and 48 agrees to 1.3e-13 V, and with the first integral to the contact tail;
    - bad inputs are refused.
- **Benchmark** (`src/verify/benchmark/tvb07.cpp`).
  - Silicon on [−1, 1] µm with N_A = 1e16 (x < 0), N_D = 2e16 (x > 0) and ohmic contacts at both ends, at equilibrium (Boltzmann, no BGN).
  - The N nodes are uniform with N even, so the step sits on the middle edge's dual face.
  - E_max is the middle edge's difference quotient, against |ψ′(0)|.
  - Inputs q, ε_si, V_T and n_i are NiTCAD's material data (I1).
  - Every run cross-checks the collocation (degree 32) against the first integral on the 64-node grid. The two differ by 1.539e-10 V, which is the predicted truncation V_T t at the contacts (1.539e-10 V). So the reference's uncertainty is measured, not assumed.
  - U_ref (truncation, relative to V_bi) is 2.1e-10, against the TB-Q4 limit of 1e-7.
  - Registered as TVB-07, tier fast. Fast runs N = 64 and 128 with V_bi only (2 comparisons). Full runs N = 64…512 (8 comparisons).
- **Results (Release, full tier):**

  | N | V_bi | ψ L∞ / V_bi | E_max relative | U_iter |
  |---|---|---|---|---|
  | 64 | 0.7339546843440443 | 4.19e-3 | 2.86e-3 | 1.7e-16 |
  | 128 | same | 1.06e-3 | 7.04e-4 | 1.4e-16 |
  | 256 | same | 2.67e-4 | 1.74e-4 | 1.4e-16 |
  | 512 | same | 6.68e-5 | 4.35e-5 | 2.1e-16 |

  - V_bi is 1 ulp from the reference 0.7339546843440442 at every level: pass (D).
  - ψ L∞ orders are 1.99, 1.98 and 1.96 (spread 0.020), so p̂ = 1.96: pass.
  - E_max orders are 2.000, 2.000 and 2.001 (spread 0.0009): pass.
  - On the finest level, ψ L∞ is 6.7e-5 ≤ 1e-4 (pass) and E_max is 4.3e-5 ≤ 1e-3 (pass). 8/8 pass, in 0.41 s.
  - Reported: the depletion edges are ε E_max/(qN), x_p = 0.2425 µm and x_n = 0.1213 µm, and ψ_j = 0.1227 V.
- **Tests.**
  - `code_verification_test.cpp` gets the fast case (2 comparisons, pass).
  - `code_verification_full_test.cpp` gets the full case (8, pass) as `verify_full.tvb07` (label `full`, 0.4 s). `verify_full.tvb04` now runs only its own tag.
  - `registry_test` expects 8 ids (was 7).
  - The verify ctests pass in Release (42/42) and Debug (42/42).
- **Time (VR10), for the owner (CV-Q19).** Serial Debug verify takes 62.6–64.5 s over two runs, against the 60 s budget.
  - On the same day the tree without step 6 measured 60.45 s, so it had no headroom left. Run-to-run spread is about ±1 s.
  - Step 6 adds 0.8 s of L0 tests (after cutting the ODE steps and collocation degrees) and 0.9 s for TVB-07's fast case.
  - TVB-09's fast case (14.7) and TVB-14's are still to come.


*Step 7, TVB-09 (MOS capacitor) and the shared first integral.*
- **Refactor, before TVB-09:** the one-sided first integral (G, the neutral bulk, Φ_b and the profile inversion) moved from `pn_junction.cpp` to `verify/reference/first_integral.{hpp,cpp}`. TVB-07 and TVB-09 both use it. The TVB-07 results and its L0 tests are bit-identical after the move (the same I3 numbers).
- **Reference** (`verify/reference/mos_capacitor.{hpp,cpp}`):
  - Kingston–Neustadter: Q_s = −sgn(t) √(2qεV_T Φ_b(u_s)), and V_G = V_FB + V_T t − Q_s/C_ox, with u_s(V_G) by Brent in double-double.
  - C_LF = C_s C_ox/(C_s + C_ox), with C_s taken analytically (its flat-band limit at t = 0). As the second method, it is also taken by complex step on a complex-double Q_s, where |t| ≥ ½. Every run checks the two against each other to 1e-10 (measured agreement is better).
  - V_FB = φ_m − χ − (E_g/2 + (V_T/2) ln(N_c/N_v)) + V_T u_b − qQ_f/C_ox.
  - L0 tests (`tests/verify/reference/mos_capacitor_test.cpp`, 6 cases, 82 assertions):
    - at flat band, Q_s = 0 and C_LF = 1/(1/C_ox + L_D/ε) to 1e-12;
    - the V_G round trip holds to 1e-25;
    - the depletion charge −√(2qεN_A V_T(t − 1 + e^−t)) holds to 1e-6. Two expectations in the test were mine and wrong, and were fixed before the commit: the first omitted the e^−t term (2.5 % at t ≈ 3), and the second took V_FB + 1.5 V for strong inversion when it is moderate (0.79 C_ox);
    - C_LF agrees with a double-double central difference of the reference's own Q_G to 1e-9, and with the complex step to 1e-10. It exceeds 0.9 C_ox in strong accumulation (V_FB − 1) and strong inversion (V_FB + 2);
    - V_FB is exact on data where it is exact by construction (u_b = −1);
    - the flat profile holds, a negative distance is refused, and bad inputs are refused.
- **Benchmark** (`src/verify/benchmark/tvb09.cpp`), as its header states:
  - p-Si with N_A = 1e17 and depth 0.5 µm; a metal gate (4.7 eV) on a lumped 10 nm SiO₂ with Q_f = 5e10 cm⁻²; an ohmic substrate.
  - The solve is Equations::equilibrium_poisson (the quasi-static state; drift-diffusion cannot reach inversion on a MOS-C, bias.hpp). It uses Boltzmann statistics and no BGN.
  - **Mesh, a design choice the criteria do not fix:** the N of the frozen ladder (64…512) are graded exponentially toward the gate, x = d(e^{aξ} − 1)/(e^a − 1) with a = 6.2. The first spacing is 0.1 nm at N = 64. The reason: the accumulation and inversion layers are under 1 nm, so a uniform mesh could not be in the asymptotic range there. h is the ξ spacing.
  - NiTCAD's side:
    - Q_s = −(Q_G + qQ_f);
    - C_LF is a Richardson-extrapolated central difference (4 and 2 mV) of NiTCAD's quasi-static Q_G at the finest level. Its uncertainty |D(h/2) − D(h)|/3 must be at most 1e-4 relative;
    - V_FB = V_G − ψ_s − Q_G/C_ox + ψ_b − qQ_f/C_ox at V_G = V_FB(ref), at every level.
  - The V_G grid is V_FB + {−1, −0.5, −0.2, −0.05, 0, 0.05, 0.2, 0.5, 1, 1.5, 2} V, swept from flat band downward and upward. Points with |Q_s| < 1e-3 of the grid's range are `not_applicable`, with the predicate named.
  - U_ref is the infinite-depth tail at x = d relative to the band bending: 1.3e-15, against 1e-6.
  - Fast tier: V_FB + 1.5 V with all four levels (7 comparisons). Full tier: the grid (37).
- **Results (Release, full tier, 0.12 s):**
  - V_FB = −0.3544794900599644 V at every level, equal to the reference: pass (D).
  - Q_s orders at every applicable V_G are 1.997–2.000, with spreads ≤ 0.002: pass.
  - Finest Q_s is 1.8e-6 to 5.6e-6, and finest C_LF is 1.5e-6 to 6.1e-6, against 1e-3: pass. The C_LF difference uncertainty is ≤ 7.9e-6.
  - At V_FB itself (Q_s = 0), the three comparisons are `not_applicable`.
  - 34 pass and 3 not applicable; overall pass.
- **Tests.**
  - The fast case is in `code_verification_test.cpp` (7, pass).
  - `verify_full.tvb09` (label `full`) freezes the 37 outcomes.
  - `registry_test` expects 9 ids.
  - The verify ctests pass in Release and Debug (43/43). Serial Debug takes 61.7–62.4 s against the 75 s budget of CV-Q19.

*Step 8, TVB-14 (full drift-diffusion manufactured solution).*
- **Construction** (`src/verify/benchmark/tvb14.cpp`, as its header states):
  - n-Si with N_D = 1e16 on [0, L]^D, L = 1 µm. Ohmic contacts at x = 0 (0 V) and x = L (V_R = 0.1 V); the sides are natural.
  - The exact state is ψ − ψ_c = V_R u + A b(u)H, φ_n = V_R u + a_n bH and φ_p = V_R u + a_p bH, where:
    - b = 4u(1 − u);
    - H = h(v)h(w), with h = 1 + (3v² − 2v³)/2 (h′ = 0 on the sides);
    - A = 0.1 V, a_n = 0.08 V and a_p = 0.15 V.
  - n and p are Boltzmann densities about NiTCAD's contact values.
  - Sources: the hook's S, G_n and G_p from the exact state in double-double. G_n = R + μ_n(∇n·∇φ_n + nΔφ_n) and G_p = R − μ_p(∇p·∇φ_p + pΔφ_p).
    - R is SRH with the material's Scharfetter lifetimes at N_D, and n_ie² = n_c p_c.
    - The mobilities are constant (μ_max).
    - Auger, radiative, BGN and Fermi–Dirac are off.
  - Errors: NiTCAD's φ_n and φ_p come from its ψ, n and p.
    - Normalisers (14.7): for ψ, max|ψ − ψ_eq| + V_T. For φ_n and φ_p, the same norm of the exact field. H1 is relative to the exact seminorm (CV-Q17).
  - Second method, run on every call: complex-step derivatives check the fields' gradients and Laplacians and ∇n (1D, 2D and 3D forms) to 1e-12.
  - Meshes: tensor N³ uniform per axis, h = L/(N − 1). The generated families use one region (the 3D generator refuses region seeds without interfaces), with h = (L^D/N)^(1/D).
- **Fast-tier family.** 14.4's families column says "1D and 2D tensor" for fast, but 14.7's frozen ladder column says "as TVB-04" (1D only in fast). 14.7 is the frozen one, so fast is 1D N = 32…256 (9 comparisons). The difference between the two sections is recorded here.
- **The sign convention** of G_n and G_p against the hook (+V G_n/R0 in the electron row, −V G_p/R0 in the hole row) is confirmed by the result: every 1D order is 1.998–2.000.
- **Results (Release, full tier, 15.6 min, 6 PARDISO threads).** 45 comparisons: 30 pass, 12 fail, 3 not evaluated.
  - Tensor orders, every field and norm: 1D 1.998–2.000, 2D 1.985–2.000, 3D 1.930–2.004. Spreads are ≤ 0.051. All 18 pass.
  - Finest L2:
    - 1D: 1.9e-5, 3.2e-5 and 1.9e-5 (ψ, φ_n, φ_p): pass.
    - 2D at 128²: 1.13e-4, 2.11e-4 and 1.18e-4: **fail**.
    - 3D at 64³: 6.7e-4, 1.43e-3 and 7.8e-4: **fail**.
    - This is the pattern of TVB-04's finding (iv). At p̂ = 2, 1e-4 needs about 2D 180² and 3D 170³ on this solution.
  - 2D generated (h_max L/128…L/16):
    - H1 orders, finest pair first: ψ 1.09, 1.27, 1.71; φ_n 1.22, 1.53, 1.90; φ_p 1.07, 1.28, 1.57. The spreads are 0.29–0.44 > 0.1, so all three are `not_evaluated` (outside the asymptotic range).
    - Finest H1 is 5.8e-4 to 1.05e-3: pass.
    - The L2 orders are 2.03–2.15 (reported).
  - 3D generated and swept (L/32…L/8, 3 levels, judged on the finest pair):
    - H1 p̂ for ψ, φ_n and φ_p: generated 1.86, 1.93, 1.52; swept 1.91, 1.97, 1.58. All six are above [0.8, 1.3]: **fail**.
    - Finest H1 is 3.4e-3 to 7.7e-3: pass.
    - The L2 orders are 2.07–2.13 (reported).
  - U_iter ≤ 6e-14 at every level, so the CV-Q8 rule is met everywhere.
  - The SRH weight max|R|/max|G_p| is 2e-6 to 6e-6 at every level (see CV-Q22).
- **Tests.**
  - The fast case is in `code_verification_test.cpp` (9, pass).
  - `verify_full.tvb14` (label `full`, TIMEOUT 3600, RUN_SERIAL) freezes the 45 measured outcomes. Since the 14.11 rulings, those are 30 pass, 12 `finding` (six (iv), six (v)) and 3 `not_evaluated`, and every result carries the CV-Q22 limitation note. The re-run after the rulings passed in 999 s (48 assertions).
  - `registry_test` expects 15 known findings (was 3).
  - **Debug time after the rulings, for the owner.** Four serial runs took 74.6, 71.2, 80.1 and 77.1 s, against 61.2–61.7 s an hour earlier on the same code apart from the finding table. Two of these exceed the 75 s budget. The increase is in TVR-1 tests that this step does not touch: `fermi_dirac` 25.6 s (was 21–22) and C-transient 22.5 s (was 15). `code_verification`, which holds every TVR-2 fast case, is unchanged at 2.4 s. The machine was under load after the full-tier runs, and the hybrid CPU schedules unpinned tests unevenly (Unit 41). Not re-run pinned; no budget change made.
  - `registry_test` expects 10 ids.
  - The verify ctests pass in Release (47/47 with the `full` label: `verify_full.tvb14` 959 s, `verify_full.tvb04` 1413 s, unchanged outcomes) and in Debug (43/43). Serial Debug takes 61.2–61.7 s against the 75 s budget of CV-Q19 (`code_verification`, which holds every TVR-2 fast case, takes 2.25 s).

**14.9 Owner rulings after the TVB-04 run (2026-11-05).**

| # | Question | Ruling |
|---|---|---|
| CV-Q12 | 2D tensor's finest L2 is 2.88e-4 against the frozen 1e-4, with both orders passing. The frozen ladders put 2D's finest h at L/63 and 3D's at L/31, against 1D's L/127, so at p̂ = 2 the frozen tolerance is reachable only on 1D. Known finding, a 256² level added after the run, or fail? | **Known finding (iv).** The criterion and ladders stay as frozen; "psi L2 finest, 2D tensor" and "… 3D tensor" are counted as `finding` with this ruling cited. The orders still gate. Nothing is widened. |
| CV-Q13 | 3D tensor (8³…64³): full, as the 14.7 ladder column reads, or nightly, as TVD-13 reads? | **Full.** |
| CV-Q14 | For studies against an exact solution, the asymptotic-range check is the order spread (≤ 0.1) alone (14.8 step 2). Confirm? | **Confirmed.** |
| CV-Q15 | Eigen SparseLU cannot hold 64³ (above). PARDISO, GMRES-ILUT, or keep Eigen? | **PARDISO for the 3D levels.** |
| CV-Q16 | With PARDISO, 64³ takes 709 s on 6 threads at about 20 GB, more than a CI runner (16 GB, 4 vCPU) holds. Full-tier ctest local only, move 64³ to nightly, or keep it in CI? | **Option 1: local only.** The full tier stays as ruled; its ctest is opt-in and not run in CI. |
| CV-Q17 | The unstructured families' H1 error: normalised by the exact ψ's discrete H1 seminorm (recommended; the potential scale A + V_T is a voltage, not a gradient norm), or literally by A + V_T? | **Recommended: the exact ψ's H1 seminorm** (owner, "choose recommended", 2026-11-05). |
| CV-Q18 | Close margins after step 5 (3D generated p̂ 1.29 against 1.3; 2D spread 0.099 against 0.1; the coarse 3D pairs above the band). Add a fourth 3D level, or keep the frozen ladders (recommended)? | **Recommended: keep the frozen ladders**; nothing added or widened (TVM2). |

**14.10 Owner ruling after step 6 (2026-11-05).**

| # | Question | Options | Ruling |
|---|---|---|---|
| CV-Q19 | Serial Debug verify is 62.6–64.5 s against VR10's 60 s. The tree before step 6 was already at 60.45 s, and TVB-09/14 fast cases are still to come. | (A) Raise VR10 to 75 s for serial Debug, stated with the measured spread. (B) Move TVB-07's fast case (and later TVB-09's) to `full`, which departs from 14.7's ladder column. (C) Keep 60 s and shorten an existing fast test (the largest are TVB-02's `fermi_dirac`, 21–22 s, and C-transient, 15 s), which touches TVR-1's frozen tests. | **A** (owner, "choose recommended", 2026-11-05). From step 6 on, the serial Debug verify budget is 75 s, measured as before (`ctest -R verify`, Debug, serial, the `full` label excluded). VR10's 60 s stays the record of TVR-1's acceptance (11.10) and is not edited. The reason: the budget had no headroom with or without TVR-2, and A changes no frozen criterion or ladder. |

**14.11 Owner rulings after step 8 (2026-11-05; owner, "choose recommended").** Nothing was widened or re-run to pass (TVM2). Each case was a miss, recorded as `fail` or `not_evaluated` until ruled; all three are ruled A.

| # | Question | Options | Ruling |
|---|---|---|---|
| CV-Q20 | TVB-14's 2D and 3D tensor finest L2 (six comparisons) miss 1e-4: 1.1e-4–2.1e-4 in 2D and 6.7e-4–1.4e-3 in 3D. Every order passes at 1.93–2.00. This is CV-Q12's situation on a second benchmark. | (A) Extend known finding (iv) to TVB-14's six comparisons. The criterion and ladder stay; the orders still gate. (B) Add finer levels: 256² is feasible, but 3D would need about 170³, which does not fit 32 GB. (C) Keep them as fail. | **A.** The six count as `finding` (iv) with this ruling cited (`registry.cpp`, `finding_iv_tvb14`), for the reason CV-Q12 gave. The orders still gate. |
| CV-Q21 | TVB-14's unstructured H1 orders: 3D generated and swept p̂ is 1.52–1.97, above [0.8, 1.3] (six fails). 2D generated is outside the asymptotic range (three `not_evaluated`), with orders falling 1.71 → 1.27 → 1.09 toward 1 as h shrinks. The 3D ladders end at L/32, where 2D's pairs were still 1.5–1.9. So this reads as pre-asymptotic super-convergence of the H1 error on these meshes, not an inconsistency. Every finest H1 passes, and every L2 order is ≈ 2. | (A) Record a known finding (v) for the six 3D fails: an H1 order above the band on the frozen 3-level 3D ladders is pre-asymptotic, and the finest-level accuracy still gates. 2D generated stays `not_evaluated`, as 14.7 prescribes. (B) Add an L/64 3D level: about 0.9M nodes and 2.7M unknowns, which does not fit 32 GB with PARDISO. (C) Keep them as fail. | **A.** Known finding (v) (`finding_v`): the six 3D H1 orders count as `finding`, and 2D generated stays `not_evaluated`. An order above the theoretical one is not evidence of a defect: the band's upper edge exists to catch a mis-specified p_th, and the tensor families gate the scheme's order cleanly. |
| CV-Q22 | SRH is on, as 14.4 requires, but on a 1 µm device with Scharfetter lifetimes (diffusion length about 50 µm) R is only 2e-6–6e-6 of the hole source. TVB-14 therefore exercises SRH but cannot detect an error in it below about that weight. | (A) Record it as a limitation of TVB-14's construction. A device near the diffusion length, or a short-lifetime material, would gate SRH, and that belongs to a later TVR. (B) Add now a second, long-device TVB-14 variant (new criteria to freeze). | **A.** Recorded as a limitation in `tvb14.cpp`'s header and as a note on every TVB-14 result. It changes no frozen item. |

**14.12 TVR-2 close-out (2026-11-06).**

| Benchmark | Fast tier | Full tier (local, MKL, 6 threads) |
|---|---|---|
| TVB-11 resistor (D) | 6 pass | — |
| TVB-12 SG exactness (D) | 8 pass | — |
| TVB-04 Poisson MMS | 3 pass | 15: 12 pass, 2 finding (iv) (CV-Q12), 1 not evaluated (3D L∞ spread, CV-Q14); 1413 s |
| TVB-07 pn junction | 2 pass | 8 pass; 0.4 s |
| TVB-09 MOS-C | 7 pass | 37: 34 pass, 3 not applicable (V_G = V_FB); 0.1 s |
| TVB-14 drift-diffusion MMS | 9 pass | 45: 30 pass, 12 finding ((iv) CV-Q20, (v) CV-Q21), 3 not evaluated (2D generated spread); 959–999 s |

- **Fast tier:** `nitcad verify --tier fast` reports 10 benchmarks: 9 pass and 1 finding (TVB-02 (ii), TVR-1's). There are no fails or run errors.
- **Full tier:** the four `full` ctests pass against their frozen outcomes (Release, 2026-11-05).
- **Tests:**
  - All 47 verify ctests pass in Release, including `full`.
  - All 43 pass in Debug.
  - Serial Debug on an idle machine takes 61.3–62.7 s (three runs), within CV-Q19's 75 s. The 71–80 s of step 8 was machine load.
- **Scope (14.1) against what was built:**
  - (a) `verify/study` builds ladders, applies the iterative-error rule and runs the order studies.
    - "GCI driving" is not used: every TVR-2 benchmark has an exact solution, so the studies judge errors against it.
    - TVR-1's GCI (`verify/reference/metrics`) applies when no exact solution exists. No TVR-2 study is of that kind.
    - The asymptotic check against an exact solution is the order spread alone (step 2, CV-Q14).
  - (b) The TVQ3 hook is the only change to existing code (step 1).
  - (c) The references have their L0 tests.
  - (d) The benchmarks are registered in the fast tier, with full-tier tests.
  - TVB-05 and 06 are deferred to TVR-4 (CV-Q3 B).
- **Recorded for later phases:**
  - The SRH weight in TVB-14 (CV-Q22): a device near the diffusion length would gate SRH.
  - The full tier is local only (CV-Q16).
  - 14.4's fast-tier families column differs from 14.7's frozen ladder column for TVB-14; 14.7 governs (step 8).

#### 15. TVR-3 (Unit 55c) proposal: the model-implementation catalogue (proposed 2026-11-06; on `verify/model-implementation`; rulings before code)

**Owner request (2026-11-06):** "start TVR-3 (55c)". TVR-3 is §9's third phase: TVB-31, an M-class check (§1.2, §4.4)
of every model in NiTCAD. This section proposes how. It also reports, before any code, places where the code's cited
source and the code may not agree (15.4). The criteria are frozen in a separate commit after these rulings and before
any code (TVM2).

**Branch.** `verify/model-implementation`, from `verify/code-verification` at `a4ed488` (PR #78, not merged): TVR-3
uses TVR-1's registry and contract and TVR-2's registry entries. It is rebased on `main` once #78 merges. Local only.

**15.1 Scope.** In:
- (a) `verify/reference/models`: a separate re-coding (I1) of each model formula in 15.3, written from the source, not
  from NiTCAD's code or structs. Base-only, as all of `verify/reference` (TVD-19; the header-boundary test applies).
- (b) The parameter register: every model parameter NiTCAD ships, typed as a decimal string from its source, with
  the source, equation or table number and the provenance grade (§3.3, PA–PD). It is data in `verify/reference`, so
  NiTCAD's material tables never feed the reference (§2.2).
- (c) The double-double functions the formulas need beyond TC1/TC2: tanh, erf, asin, log1p and ln(1 + eˣ), and the
  dilogarithm (for F₁). Each gets L0 tests against DLMF identities and published values, as TC1–TC9 did.
- (d) TVB-31's entries in the registry, one per model (TVB-31.01 …), in the `fast` tier, with their tests.

Out:
- The Fermi–Dirac integrals F_{1/2}, F_{3/2} and their inverse (TVB-02), the Bernoulli function (TVB-03) and the
  constants (TVB-01). These are already verified in TVR-1.
- Physics as assembled: the SG flux, the heterojunction and Schottky boundary terms, the BTBT path search, the
  gate-tunnelling coupling. These are TVR-4/TVR-5 benchmarks.
- Any change to NiTCAD's code, tables, tests, gates or tolerances. A miss is recorded as a finding, and any fix is
  proposed separately (as Unit 56 did for findings (i)–(iii)).

**15.2 What an M check compares.** For each model and each point of its grid (15.5):
- **The value** NiTCAD's physics function returns, against the reference evaluated in double-double and rounded.
- **Every partial derivative NiTCAD returns** (`d_dE`, `d_dn`, the thermal slopes, …), against the reference's
  analytic derivative. Each reference derivative is checked at L0 by complex step (TC8) on the reference's own formula,
  so a Newton-facing derivative is verified against the source and not only against NiTCAD's FD-Jacobian gates (which
  are I0, K1).
- **The parameters** (the register, (b)), bit-equal to the correctly rounded parse of the source's printed decimal.
- **The source's printed values**, where it prints computed results (e.g. Couderc's n_i at stated temperatures), to
  half a unit in the last printed digit (the NPARC ruling of 11.5).

The benchmark side calls NiTCAD's public physics functions with parameters from the register, so that a formula
error and a table error are told apart:
- formula check: NiTCAD's function with the register's parameters vs. the reference;
- table check: NiTCAD's shipped material set vs. the register.

**15.3 The models (the catalogue's "every model").** Grades are what DECISIONS.md and the headers record today; step 1
re-reads each source. "Legacy" means the value or form came from the legacy code and no source for it has been read in
this repository. That is weaker than PC, and is listed as a grade of its own (MI-Q5).

| Entry | Model (function) | Source the code cites | Class | Grade today | Grid (stated range) |
|---|---|---|---|---|---|
| 31.01 | Varshni E_g(T) (`band_gap_eV`, legacy set) | Varshni form; silicon values legacy | M | legacy | T 1–500 K (header) |
| 31.02 | Pässler E_g(T) (default set) | Couderc et al. 2014 eq (10), Table I | M | PA (Unit 35) | 1–400 K |
| 31.03 | DOS mass fits, N_c(T), N_v(T) | Couderc eqs (4), (5), (18), (19) | M | PA | 1–400 K |
| 31.04 | n_i(T), E_c − E_i, d ln n_i/dT, dE_g/dT, d ln N/dT (`thermal.hpp`) | derived (PD) | E (formula) | PD | 1–500 K |
| 31.05 | Slotboom BGN and n_ie | Slotboom (de Graaff); values legacy | M | legacy | N 1e14–1e21 cm⁻³ |
| 31.06 | Caughey–Thomas μ(N, T) and dμ/dT | Caughey–Thomas 1967; values legacy | M | legacy | N 0, 1e13–1e21; T 200–500 K |
| 31.07 | Masetti μ(N, T), dμ/dT (B, P) | Masetti 1983 eq (1), (4), Table I | M | PB (B), PC (P μ_max) | N 0, 1e13–1e21 |
| 31.08 | Canali μ(E), dμ/dE | Canali 1975; v_sat, β legacy | M | legacy | E 0, 1e-2–1e6 V/cm |
| 31.09 | Lombardi μ and its four partials | Lombardi 1988 | M | PB | E⊥ 1–1e7 V/cm; N; T |
| 31.10 | Carrier-temperature mobility, E_eff (Unit 52) | derived from Canali (PD) | E (formula) | PD | a −1e-3…1e6; β ∈ {1, 2, 1.5} |
| 31.11 | Scharfetter lifetimes τ(N) | legacy | M | legacy | N 0, 1e13–1e21 |
| 31.12 | SRH R and partials (Boltzmann and FD E) | Shockley–Read, Hall 1952 | M (form) | standard form; not read | n, p over 1e-10…1e21 with np/E across 1 |
| 31.13 | Auger R and partials | Dziewior–Schmid values (legacy) | M | legacy | as 31.12 |
| 31.14 | Radiative R (B per material) | Ioffe; nextnano (GaAs) | M | PC | as 31.12 |
| 31.15 | Interface-trap occupancy and rate | Shockley–Read form (Unit 15b) | M (form) | PD from SRH | x, n, p, n₁, p₁ |
| 31.16 | Trap-band quadrature (`trap_band_levels`) | 6-point Gauss–Legendre, panels ≤ kT | E (in-model numerics) | PD | bands up to the gap, 200–500 K |
| 31.17 | Incomplete ionization N_D⁺, N_A⁻ | standard form; 45 meV, g = 2, 4 legacy | M | legacy | η −40…40, depth 0–0.3 eV |
| 31.18 | Altermatt–Schenk E_dop(N), b(N) | Schenk et al. SISPAD 2006 eqs (1), (5), Table I | M | PA (SISPAD), JAP not read | N 1e14–1e21 |
| 31.19 | Neutral equilibria (Boltzmann, FD, ionized) | root of a closed form | E (in-model numerics) | PD | net doping ±1e21; compensated |
| 31.20 | vOdM α(E), dα/dE; γ(T), dγ/dT | vOdM 1970; γ "Sentaurus' form" | M | legacy (values), unread (γ) | E 1e4–1e7 V/cm and branch switch; T 200–500 K |
| 31.21 | Kane G(F), dG/dF | Hurkx–Klaassen–Knuvers 1992 Table I (values) | M | legacy (attribution unread) | F 1e4–1e7 V/cm |
| 31.22 | Direct-gap WKB BTBT: κ(δ), its antiderivatives, segment integrals, path rate, eq (8) | Esseni et al. 2017 §2.1 eqs (8), (9), (11), (12) | M | not recorded | u > 1; δ ∈ [0, 1]; F |
| 31.23 | WKB segment action, transmission | Gehring 2003 §3.5.1 (closed-form integral) | E (formula) | PA (open thesis) | a, b across 0 and a = b |
| 31.24 | Tsu–Esaki supply N(E) and partials | Gehring §3.2, 3.3.1 | M (form) | PA | (E − F)/kT −60…60 |
| 31.25 | Tsu–Esaki current (4-point GL quadrature) | Gehring | E (in-model numerics) | PA | MI-Q8 |
| 31.26 | Thermionic emission velocity (both forms), dv/dT | Richardson law, mass from N (PD) | E (formula) | PD | T 200–500 K |
| 31.27 | MLDA λ, f(z), mean factor and its series | Paasch–Übensee; Hänsch (via Ieong, Entner) | M; mean factor E (formula) | PB (form), PD (mean) | z/λ 0–10; a/λ 1e-8–1e3 |
| 31.28 | κ(T) power law, ρc | legacy; Ioffe; NIST-JANAF | M | legacy / PC | 200–600 K |
| 31.29 | Thermopower factor h(x), F₋₁, F₀, F₁ | Wachutka 1990; Lundstrom §3.4; F_j exact | E (special functions) | PD + DLMF | x −750…750 |
| 31.30 | SiO₂ data: ε_r, gap, offset, masses, κ | legacy 3.9; Gehring Table 5.1; owner ruling (masses) | M (table) | legacy / PA / choice | — |
| 31.31 | Material sets: Ge, GaAs, In₀.₅₃Ga₀.₄₇As, 4H-SiC, Al_xGa₁₋ₓAs(x) | Adachi; Ioffe; legacy | M (table) | legacy / PC | x ∈ [0, 0.45] |

Classes in this table: M is §1.2's empirical model. "E (formula)" is an exact consequence derived here (PD): by the
contract of 11.4 it needs a second method, so each such reference has two:
- the closed form against quadrature (TC3/TC4) of the integral it solves (31.23, 31.27 mean, 31.29), or
- against the root (TC5) of the equation it inverts (31.10, 31.19), or
- against complex-step differentiation of the reference's value (the derivatives in 31.04, 31.26).

"E (in-model numerics)" is a numerical approximation inside a model (31.16, 31.19, 31.25), judged against a
double-double reference integral or root.

**15.4 Found before code (from the headers and from recollection; each to be confirmed by reading the source in step
1, so none is yet a finding).**
- **(F1) Slotboom offset.** `bandgap_narrowing.hpp` adds −√½ to Slotboom's form so that ΔE_g is 0 at N₀. The
  published form (ΔE_g = E₀[ln(N/N₀) + √(ln²(N/N₀) + ½)]) has no such term. Above N₀ the difference is a constant
  E₀√½ = 4.89 meV, about 9.9 % in n_ie at 300 K. The header records it as the legacy's, so it is a documented
  departure, not a defect. Separately, E₀ = 6.92 meV and N₀ = 1.3e17 cm⁻³ are, as I recall, Klaassen, Slotboom and
  de Graaff's 1992 refit (Solid-State Electron. 35, 125), not the 1976–77 Slotboom–de Graaff values (9 meV, 1e17).
- **(F2) Canali β.** `field_mobility.hpp` attributes silicon's β = 2 (electrons) and 1 (holes) to Canali et al. 1975.
  As I recall, Canali's fit is β_n = 1.109 (T/300)^0.66 and β_p = 1.213 (T/300)^0.17, and β = 2 and 1 are Caughey and
  Thomas's 1967 field-dependent fit. Canali's v_sat law, as I recall it (1.53e9 T^−0.87 and 1.62e8 T^−0.52 cm/s),
  gives 1.0704e7 cm/s for electrons at 300 K, which agrees with the shipped 1.07e7. It gives 8.345e6 cm/s for holes,
  against the shipped 8.37e6 (0.3 %).
- **(F3) Caughey–Thomas silicon values.** The electron set (92, 1360, 1.3e17, 0.91) is, as I recall, Baccarani and
  Ostoja's (1975). The hole set (47.7, 495, 6.3e16, 0.76) is Caughey and Thomas's. The temperature exponents (−2.33,
  −2.23) are Arora, Hauser and Roulston's (1982). None of these sources is cited or has been read.
- **(F4) The impact-ionization temperature factor** cites "Sentaurus' vanOverstraetendeMan model". That manual is
  proprietary and not citable as provenance. Its primary source needs finding, or the factor is graded a NiTCAD choice.
- **(F5) The Kane pair** (3.5e21 cm⁻³ s⁻¹, 1.03e8 V/cm) is attributed to Hurkx, Klaassen and Knuvers 1992, Table I.
  Its form and units have not been checked against that table.
- **(F6) Masetti's** μ_max(T) exponents are stated in the header not to be Masetti's. They are a NiTCAD choice, and
  phosphorus μ_max = 1414 is PC (C-10). Both are already recorded.

(F1) and (F6) are documented departures. (F2)–(F5) are attribution questions. Under MI-Q4, a documented departure is
checked twice:
- the code against its own documented form (an implementation check, which should pass);
- the code against the published form (a miss here becomes a finding with its attribution).

**15.5 Grids and the criterion.**
- **Grids.** Each model's grid spans the source's stated range (table 15.3), log-spaced where the variable is, at
  ≥ 64 points per decade. Each grid also includes:
  - every switch point of the code: vOdM's branch field, Slotboom's N₀, MLDA's series switch, the F_j switch at
    |x| = 1, the regularization a₀ of 31.10, and each underflow cut;
  - the double neighbours on each side of every switch point (`nextafter`);
  - zeros of the value where the model has them (np = E for 31.12–31.15).
  Outside the stated range, points are evaluated and reported as `not_applicable` ("outside the source's range",
  §4.4 M), never as passes.
- **Criterion (MI-Q2).** §4.4 fixes M at 1e-13 relative. That cannot hold where the formula itself is
  ill-conditioned. For example, exp(−B/E) at B/E = 700 amplifies a 1-ulp input error 700 times (1.6e-13), and np − E
  near equilibrium cancels. TC1's ruling scaled the exp identities by their condition number for the same reason. The
  proposal:

      |q_sim − q_ref| ≤ max(1e-13 |q_ref|, 64 ε Σᵢ |xᵢ ∂q/∂xᵢ|)

  - The sum runs over the inputs, so it is the absolute form of the condition number κ|q|, from the reference's own
    derivatives.
  - It is 1e-13 relative wherever κ ≤ 7, and it holds at a zero of q (where relative error is undefined).
  - The same rule applies to each partial derivative. Its κ uses the reference's second derivatives, taken by
    double-double central differences of the analytic first derivatives (a tolerance scale only, not a judged
    quantity).
- **Documented cuts** (MI-Q7). Examples: α = 0 where γB/E > 700, G = 0 where B/F > 700, ΔE_g = 0 at N ≤ N₀. These are
  judged as cuts. NiTCAD's 0 passes when the reference's value at that point is below the documented threshold's
  value. The miss of a cut that sits in the wrong place is caught by the switch-point neighbours.
- **In-model numerics** (31.16, 31.19, 31.25). The criteria are in MI-Q8.

**15.6 Order of work (one commit each, L0 tests before the benchmark that uses them; mutation runs per commit as in
TVR-1/2).**
1. Sources and the register. Read each source in 15.3 (primaries where open, otherwise two independent open
   reproductions for PB: TU Wien theses, Xyce, Genius, DEVSIM, Allpix Squared, Silvaco's public manual, nextnano,
   Ioffe). Write the register with grades, and record 15.4's items as confirmed or withdrawn. Documentation and data
   only.
2. Double-double tanh, erf, asin, log1p and ln(1 + eˣ), and Li₂, with L0 tests.
3. Band set and narrowing: 31.01–31.05.
4. Mobility: 31.06–31.10.
5. Recombination and traps: 31.11–31.16.
6. Statistics and ionization: 31.17–31.19.
7. Generation: 31.20–31.22.
8. Tunnelling, emission, MLDA, thermal: 31.23–31.29.
9. Material tables: 31.30, 31.31, and every shipped set against the register.
10. Close-out: outcomes, the provenance summary (R-3's register surfaced in the report), Debug re-timed, full suites.

**15.7 Owner questions (recommendations in the last column).**

| # | Question | Recommendation |
|---|---|---|
| MI-Q1 | Scope: table 15.3 as "every model", with the exclusions of 15.1 (TVB-01–03 already done; assembled physics in TVR-4/5)? | Yes |
| MI-Q2 | Criterion: amend §4.4 M from a fixed 1e-13 relative to 15.5's max(1e-13 \|q\|, 64 ε Σ\|xᵢ ∂q/∂xᵢ\|), for values and partials alike? (A) As proposed. (B) Keep a fixed 1e-13 and report the ill-conditioned points as `not_evaluated` | A. It is 1e-13 wherever the formula is well-conditioned, and B would leave the underflow edges and equilibrium unjudged |
| MI-Q3 | Judge every partial NiTCAD returns against the reference's analytic derivative, itself checked by complex step? | Yes. They are what Newton uses, and today only I0 FD gates check them |
| MI-Q4 | Where the code documents a departure from its source (F1, F6), check against both the documented form (pass expected) and the published form (a miss is a finding)? | Yes |
| MI-Q5 | Provenance grades: (A) add a grade **L** ("legacy: no source read"), treated as PC: such an entry cannot pass and is a finding until step 1 reaches PB or better; NiTCAD's own choices (owner rulings, e.g. the SiO₂ masses, EB β = 2, Masetti's T exponents) get a tag **C** and are pinned to the ruling, not verified. (B) Treat legacy values as PB | A. B would credit the legacy code as a source |
| MI-Q6 | Sources: step 1 aims at PB from open reproductions, and the owner supplies primaries where they can (TVR0-Q2, still open). Priority list: Caughey–Thomas 1967, Baccarani–Ostoja 1975, Arora et al. 1982, Canali et al. 1975, Slotboom–de Graaff 1976/77 and Klaassen et al. 1992, Dziewior–Schmid 1977, van Overstraeten–de Man 1970, Hurkx et al. 1992, Masetti et al. 1983, Lombardi et al. 1988, Paasch–Übensee 1982, Hänsch et al. 1989, Altermatt et al. 2006, Esseni et al. 2017, Varshni 1967, Green 1990. | Yes. Grades are raised whenever a primary arrives |
| MI-Q7 | Grids per 15.5: the source's range at ≥ 64 points per decade, switch points and their neighbours, zeros; out of range `not_applicable`; documented cuts judged as cuts? | Yes |
| MI-Q8 | In-model numerics: (a) trap-band quadrature (31.16) and the neutral-equilibrium roots (31.19) under 15.5's criterion against double-double references; (b) the Tsu–Esaki current quadrature (31.25) **measured and reported without pass/fail** (as §4.4 H), with a threshold set by the owner afterwards. Its integrand has kinks where a turning point crosses a segment end, so a frozen 4-point GL accuracy has no stated bound to judge against | Yes |
| MI-Q9 | Tier and cost: all of TVB-31 in `fast`. TVR-3 adds at most 8 s to the serial Debug verify (61–63 s of the 75 s of CV-Q19). Beyond that, the dense grids move to `full` and `fast` keeps every switch point and one point per decade | Yes |
| MI-Q10 | Findings: no code, table or attribution in NiTCAD is changed in TVR-3. Each confirmed item of 15.4, and each miss, goes to the owner as a finding with a fix proposal, as §12 did | Yes |
| MI-Q11 | Unit and branch: Unit 55c on `verify/model-implementation`, stacked on PR #78 and rebased on `main` after it merges. Order of work as 15.6 | Yes |

**15.8 Rulings (owner, "choose recommended", 2026-11-06).** MI-Q1 through MI-Q11 are approved as recommended:
- The scope is table 15.3, with 15.1's exclusions.
- §4.4's M criterion is amended for TVR-3 to 15.5's condition-scaled form, for values and partials alike.
- Every returned partial is judged.
- Documented departures (F1, F6) are checked against both forms.
- Grade **L** (legacy, no source read) is added and treated as PC. Tag **C** (NiTCAD's own choice) is pinned, not
  verified.
- Step 1 aims at PB from open reproductions. The primaries of MI-Q6 are welcome whenever supplied.
- The grids are 15.5's.
- Trap-band quadrature and the neutral-equilibrium roots are judged. The Tsu–Esaki quadrature is reported only.
- TVB-31 goes in `fast`, with an 8 s Debug allowance.
- No NiTCAD change is made in TVR-3.
- Unit 55c, on this branch.

**15.9 Criteria (frozen 2026-11-06, before code; TVM2: never widened after a run).** Notation: ε = 2⁻⁵²; q_ref is
the double-double reference rounded to double; S(q) = Σᵢ |xᵢ ∂q/∂xᵢ| over the model's inputs, from the reference's
derivatives.

| Check | Criterion | Outcome rules |
|---|---|---|
| Value (every entry) | \|q_sim − q_ref\| ≤ max(1e-13 \|q_ref\|, 64 ε S(q)) at every grid point inside the source's range | Any point over: fail (or finding, under grade L/PC, per 11.4's contract). Out of range: `not_applicable`, value reported |
| Partials | The same rule for each returned partial, with S of that partial from double-double central differences of the reference's analytic derivative. Each analytic derivative is first checked at L0 against complex step to 64 ε S | As for values |
| Documented cuts | NiTCAD returns exactly 0 where the code documents a cut. That passes if \|q_ref\| there is ≤ \|q_ref\| at the cut's threshold. Both `nextafter` neighbours of every cut and switch point are judged by the value rule | A non-zero inside a cut, or a zero outside one: fail |
| Parameter register | Each shipped parameter is bit-equal to the correctly rounded parse of the source's printed decimal | A mismatch is a finding, recorded with its source line. Grade L or PC: finding at best. Tag C: pinned to the ruling's value |
| Printed values | The reference, and NiTCAD, reproduce the source's printed computed values to half a unit in the last printed digit | Miss: finding |
| Documented departures (F1, F6) | The documented form: the value rule. The published form: the value rule | The published-form miss is a finding with its attribution |
| Grid density | ≥ 64 points per decade (or per unit for linear variables, ≥ 256 points over the range), plus switch points, their neighbours, and zeros | — |
| 31.16 trap band, 31.19 roots | The value rule against a double-double reference: an adaptive Gauss–Kronrod integral (TC3) of the occupancy, and a Brent root (TC5) | As for values |
| 31.25 Tsu–Esaki current | Relative difference from a double-double adaptive integral, reported per point; no pass/fail until an owner threshold | `diagnostic` only (never pass, fail or finding) |
| Cost | TVB-31 adds ≤ 8 s to the serial Debug verify (idle machine). Beyond that, MI-Q9's fast/full split applies, with no change to the criteria | — |

References carry U_ref from their double-double evaluation. By the TB-Q4 rule, a reference whose U_ref exceeds 1e-3 of
the bound at a point is `not_evaluated` there.

**15.10 Step 1: reading log (in progress, 2026-11-06).** The pages were read by automated fetch of open copies.
"Repro" means an open reproduction, not the primary. No grade is raised to PA from these reads.

| Source read | What it gives (as printed) |
|---|---|
| Tidy3D `SlotboomBandGapNarrowing` (docs.flexcompute.com) | ΔE_g = V₁[ln(N/N₂) + √(ln²(N/N₂) + C₂)] for N ≥ min_N (1e15), else 0. **No constant is subtracted.** Example silicon set: 6.92e-3 V, 1.3e17, 0.5. Cites *Unified apparent bandgap narrowing in n- and p-type silicon*, Solid-State Electron. 35, 125 (1992) |
| COMSOL 6.3 Semiconductor Module, Band Gap Narrowing | Its silicon library uses "the updated material properties due to Klassen and others … rather than the original properties". The equation is an image and was not read |
| Palankovski, PhD thesis, TU Wien 2000, Table 3.15 | Slotboom silicon: E_ref 0.009 eV, C_ref 1e17 cm⁻³; cites Slotboom and de Graaff, Solid-State Electron. 19, 857 (1976) |
| Allpix Squared 3.3 manual, mobility | Canali, silicon. Electrons: v_m = 1.43e9 T^−0.87 cm/s (Canali eq 2a); the Jacoboni–Canali table gives 1.53e9 T^−0.87; β_e = 2.57e-2 T^0.66. Holes: v_m = 1.62e8 T^−0.52; β_h = 0.46 T^0.17 (Jacoboni et al., Solid-State Electron. 20, 77 (1977), Table 5). Masetti B and P sets agree with NiTCAD's to the printed digits. Arora: μ₀ = 7.40e8 T^−2.33 (e), 1.36e8 T^−2.23 (h) |
| Allpix Squared 3.3 manual, impact ionization | vOdM coefficients from the 1970 abstract: E₀ 4.0e5 V/cm; e 7.03e5 /cm, 1.231e6 V/cm; h low 1.582e6, 2.036e6; h high 6.71e5, 1.693e6. γ(T) is the tanh ratio "following the Synopsys Sentaurus TCAD user manual", with ħω 0.063 eV. "T₀ … is never stated explicitly [in the original work]; a value of T₀ = 300 K is assumed" |
| Palankovski, Tables 3.5, 3.8, 3.20, 3.27, 3.31, 3.36 | κ = 148 W/(K m) (T/300)^−1.3 (refs Ng, Sze). Varshni Si: E_g0 1.1695 (reported 1.17), α 4.73e-4, β 636 (refs Singh, Sze). Green's N_C,V(T) polynomial (Green, J. Appl. Phys. 67, 2944 (1990)) gives 2.86e19 and 3.10e19 at 300 K (my sum). DD high-field β: Si 2.0 (n), 1.0 (p), uncited. Quay v_sat: 1.0e7 and 7.04e6 cm/s. Auger Si: 2.8e-31, 9.9e-32 cm⁶/s, uncited |
| Gritsch, PhD thesis, TU Wien, §4.7 | The "Scharfetter relation" τ = τ_min + (τ_max − τ_min)/(1 + (N/N_ref)^γ), cited to Fossum and Lee (1982) and Fossum et al. (1983). No silicon values printed |
| Jiao et al., J. Semicond. 34, 092002 (2013), Table 1 | Si BTBT (A, B, P): direct 1.35e20 cm⁻³s⁻¹, 101 MV/cm, 2; indirect 3.29e15, 23.8 MV/cm, 2.5; Sentaurus (Hurkx) defaults 4e14, 19 MV/cm, 2.5 |
| Ayalew, PhD thesis, TU Wien, §3.3.3 | The Canali form with a T-law; SiC numbers only. Cites Canali et al., IEEE TED 22, 1045 (1975) |

Status of 15.4's items after these reads:
- **F1 confirmed.** One reproduction prints the published form (Tidy3D), and it has no −√½ term. COMSOL's equation
  is an image, and Palankovski's was not reported. This makes it **PC** for the form, with a second reproduction still
  needed. On the values, three reproductions agree (**PB**): 6.92 meV and 1.3e17 are the 1992 refit (Tidy3D, COMSOL),
  and 9 meV and 1e17 are the 1976 values (Palankovski). NiTCAD's header attributes the 1992 values to "Slotboom (de
  Graaff)", which is an attribution finding. The legacy offset stays a documented departure (MI-Q4).
- **F2 confirmed (PB).**
  - Canali's β at 300 K is 0.0257·300^0.66 = 1.109 for electrons and 0.46·300^0.17 = 1.213 for holes (Allpix).
  - β = 2 and 1 appear in Palankovski uncited, as the generic DD high-field exponents, not as Canali's.
  - The electron v_sat of 1.07e7 is the Jacoboni–Canali law (1.53e9·300^−0.87 = 1.0704e7). Canali's own eq 2a gives
    1.43e9·300^−0.87 = 1.0005e7.
  - The hole v_sat from the shared law, 1.62e8·300^−0.52 = 8.345e6, differs from the shipped 8.37e6 by 0.30 %. No read
    source prints 8.37e6.
- **F3 open.** No source read prints NiTCAD's Caughey–Thomas sets. One search result attributes a 1360/92/1.3e17/0.91
  electron set to Baccarani and Ostoja (not read). The T exponents −2.33 and −2.23 are Arora's (Allpix), applied to
  μ_max. In Arora's model they apply to μ₀ in a different form.
- **F4 confirmed.** Every reproduction read takes γ(T) from the Sentaurus manual, and none from a primary. Allpix
  states that the original vOdM paper gives no T₀. The factor is a tag **C** (a NiTCAD choice) until a primary is
  found. vOdM's coefficients agree with Allpix's reproduction of the 1970 abstract (one reproduction, PC). NiTCAD's
  electron switch at 5e5 V/cm is immaterial, since both branches are equal.
- **F5 strengthened.** The pair 3.5e21 cm⁻³s⁻¹ with 1.03e8 V/cm appears in none of Jiao et al.'s three silicon sets.
  Sentaurus's Hurkx defaults are 4e14, 1.9e7 V/cm and P = 2.5, a different form from NiTCAD's F². The attribution
  to Hurkx et al.'s Table I is unsupported by every source read: grade **L**.
- **F6 unchanged.** The Masetti sets match Allpix (a second reproduction for phosphorus). This could raise P's μ_max =
  1414 from PC to PB: Allpix and the earlier Palankovski Table 3.24 make two.
- **Other entries:**
  - Varshni silicon: one reproduction (Palankovski, citing Sze and Singh), PC.
  - Legacy N_c and N_v are Green's 1990 polynomial at 300 K (one reproduction, PC).
  - Auger: values reproduced uncited. Dziewior–Schmid unconfirmed (L for the attribution, PC for the values).
  - κ: 1.48 W/(cm K) is reproduced, but the exponent is printed as −1.3, not 1.33 (a 0.03 difference, origin unknown,
    L).
  - Scharfetter lifetimes: the form is cited to Fossum et al. The values 1e-5, 3e-6 and 5e16 are not found (L).

**Second batch (2026-11-06).**

| Source read | What it gives (as printed) |
|---|---|
| Esseni, Pala, Palestri, Alper, Rollo, Semicond. Sci. Technol. 32, 083005 (2017), open access (CC BY 3.0), pp. 4–5; PDF rendered locally | eq (8): G_T = e²\|F\|²√m_r / (18π ħ² √E_G) · exp(−π √m_r E_G^{3/2} / (2ħe\|F\|)). eq (9): κ = √(m_r E_G (1 − α²))/ħ, α = −m₀/(2m_r) + 2√((m₀/(2m_r))((E − E_V)/E_G − ½) + m₀²/(16m_r²) + ¼). eq (11): G_T(E) = \|dE_V/dx\| e/(36ħ) (∫dx/κ)⁻¹ [1 − exp(−k_m² ∫dx/κ)] exp(−2∫κ dx) (f_c − f_v). eq (12): k_m² = min(2m_v(E_MAX − E)/ħ², 2m_c(E − E_min)/ħ²). m_r⁻¹ = m_c⁻¹ + m_v⁻¹ |
| Genius TCAD Open, `src/material/Si` (raw files) | BGN: V0 6.92e-3 eV, N0 1.3e17, CON 0.5, EgNarrow = V0(x + √(x² + CON)), **nothing subtracted**. Auger 2.8e-31, 9.9e-32 (uncited). NSRHN = NSRHP = 5e16 but TAUN0 = TAUP0 = 1e-7. NC300 2.86e19, NV300 3.10e19 (comment: Green 1990). Dopants: B 0.045 eV, g 4; P 0.045 eV, g 2. Analytic mobility 55.24/1429.23/1.072e17/0.73 (e), 49.70/479.37/1.606e17/0.70 (h); β 2 and 1, uncited. OP_PH_EN 0.063 eV. τ_w 6.8e-13 and 2.0e-13 s. Its BTBT A and B were seen only in a fetch summary (unverified) |
| Xyce reference guide tables and `N_DEV_Material*.C` | Caughey–Thomas Si as Genius's analytic set. Arora exponents −2.33 (exn2) and −2.23 (exp2). Auger 2.8e-31, 9.9e-32. Lifetimes from Weaver–Nasby (3.95e-4 s and 3.52e-5 s, N_ref 7.1e15). E_g(T) = 1.17 − 4.73e-4 T²/(T + 636), cited to Thurmond (1975). Slotboom "from the Medici manual": 9e-3, 1e17, 0.5. N_v 2.66e19 |
| Triebl, PhD thesis, TU Wien 2012, §4.3.1 | "The energy relaxation time is often assumed to be constant, with typical values τ_E,n = 0.35 ps [Vasicek 2009] and τ_E,p = 0.4 ps [Palankovski 2000]" |

Effect:
- **31.22 is PA.** NiTCAD's κ, α (with u = m₀/(2m_r)), the path rate and k_m² match eqs (9), (11) and (12)
  term for term. Its uniform-field rate matches eq (8) (`band_to_band.cpp:157`, in SI with |dE_V/dx| in J/m). One
  printed sign is inconsistent: eq (11) has (f_c − f_v), while eq (2) has (f_v − f_c) with f_v = 1, f_c = 0. NiTCAD
  takes the latter, which makes it pure generation. This is recorded as a misprint, not a finding.
- **F1's form reaches PB.** Genius prints Tidy3D's form with nothing subtracted. Xyce attributes 9 meV / 1e17 to the
  Medici manual's Slotboom.
- **Auger values reach PB** (three reproductions). None of them cites Dziewior–Schmid.
- **Legacy N_c and N_v reach PB** (Green 1990, via Genius and Palankovski). The Varshni set reaches PB (Thurmond 1975,
  via Xyce and Palankovski).
- **F3 is still open after four codes and two theses.** None prints 92/1360/1.3e17/0.91 or 47.7/495/6.3e16/0.76, so
  the grade stays L.
- **The lifetimes** 1e-5 and 3e-6 are still not found (L). N_ref = 5e16 is in Genius, uncited (PC).
- **τ_w** is PA as the code cites it (Triebl's "typical values").

**The register** is `include/NiTCAD/verify/reference/model_parameters.hpp`: 104 rows typed from these reads.
- Each row carries its grade (MI-Q5 L and C included) and its source locator.
- Its silicon, Masetti, Lombardi, Altermatt, Canali-law, vOdM, Kane, thermal, SiO₂ and MLDA rows cover every
  shipped silicon model parameter, except the thermopower exponents (−½, a theory value, 31.29) and fields shipped as 0 meaning "not modelled" (radiative B, Richardson constants).
- Fields with no source read (ε_r, χ, the Caughey–Thomas sets, lifetimes, the Kane pair, ρc, SiO₂ ε_r) are rows with
  an empty value and grade L.
- The non-silicon sets (31.31) and the radiative coefficients (31.14) are read in step 9, where they are used. This
  is a reordering of 15.6 that changes no criterion.

Remaining source acquisition (open copies exhausted; the MI-Q6 primaries would settle these): Caughey–Thomas 1967 and
Baccarani–Ostoja 1975 (F3), the lifetime set, Dziewior–Schmid 1977, the Kane pair's origin (F5), and the γ(T) primary
(F4).

Superseded list (first batch): Caughey–Thomas and Baccarani–Ostoja values; Dziewior–Schmid; the hole v_sat 8.37e6; the
Scharfetter values; the Kane pair's origin; Esseni et al. 2017; the dopant levels and degeneracies; radiative
coefficients; the SiO₂ ε_r; the non-silicon material sets; Triebl's τ_w (§3.3 provenance for 31.10's inputs). Then the
register file itself.

**15.11 Step 2: L0 criteria for the new double-double functions (frozen 2026-11-06, before code; owner "choose
recommended and proceed").**
- **Functions:** `log1p`, `log1p_exp` (ln(1 + eˣ)), `tanh`, `erf`, `erfc`, `sin` and `cos` (for |x| ≤ 1), `atan`,
  `atan2`, `asin` and `li2` (real x ≤ 1). They go in `verify/reference/special_functions`, and depend on
  `double_double` only.
- **Methods** (sources DLMF §4.24, §4.26, §7.6, §7.9, §25.12; PA, open):
  - log1p: 2 atanh(x/(2 + x)) by series for |x| < 1/8.
  - erf: the series of positive terms 2/√π e^{−x²} Σ 2ⁿx^{2n+1}/(2n+1)!! (DLMF 7.6.2) for x < 3, and erfc by the
    even-part continued fraction √π e^{x²} erfc x = 2x/(2x² + 1 − 1·2/(2x² + 5 − 3·4/(2x² + 9 − ⋯))) (DLMF 7.9.2, read
    2026-11-06) for x ≥ 3.
  - atan: Newton on tan from std::atan, with sin and cos by Taylor.
  - asin(x) = atan2(x, √((1 − x)(1 + x))).
  - li2: the series (DLMF 25.12.1), the Landen identity (25.12.3), reflection (25.12.6) and inversion (25.12.4),
    all read 2026-11-06.
- **Criteria** (u_dd = 2⁻¹⁰⁴; U_ref is the quadrature's own error estimate, added to the bound):

| Check | Bound |
|---|---|
| erf(1) against OEIS A099286's 40 digits; tanh(1) against A073744's; exact comparison (TC1's method) | 4 u_dd relative |
| li2(−1) = −π²/12, li2(½) = π²/12 − ½ ln²2, li2(1) = π²/6; atan(1) = π/4; asin(½) = π/6; sin(π/6) = ½; cos(π/4) = √2/2 | 8 u_dd relative |
| erf(x) against qag15 of (2/√π) e^{−t²} over [0, x], 64 points on (0, 6] | 16 u_dd relative + U_ref |
| erfc(x) against qag15 of (2/√π) e^{−t²} over [x, x + 15], 64 points on [3, 25] (amended before any run: the first draft said [3, 26], but erfc(26) ≈ 5e-296 is below double-double's normal range of 1e-280, `double_double.hpp`) | 16 u_dd relative + U_ref |
| atan(x) against qag15 of 1/(1 + t²) over [0, x], 64 points on (0, 20] | 16 u_dd relative + U_ref |
| li2(x) against −qag15 of ln(1 − t)/t over [0, x], 64 points on [−20, 0.99] | 16 u_dd relative + U_ref |
| log1p(x) against x − x²/2 + x³/3 − x⁴/4, for \|x\| ≤ 1e-9 | 4 u_dd relative |
| expm1(log1p(x)) = x, for x ∈ (−0.9, 10] log-spaced both signs | 16 u_dd relative |
| log1p_exp(x) − log1p_exp(−x) = x, for \|x\| ≤ 600 | 8 u_dd · max(1, \|x\|) absolute |
| tanh(2x) = 2t/(1 + t²), with t = tanh(x), for \|x\| ∈ [1e-20, 40] | 16 u_dd relative |
| sin²y + cos²y = 1 for \|y\| ≤ 1; sin(asin x) = x for \|x\| ≤ 1; atan(x) + atan(1/x) = π/2 for x > 0 | 8, 16, 8 u_dd |

A miss is brought to the owner. It is not met by widening a bound (TVM2).

**15.12 Step 2 as built (2026-11-06).** The files are `include/NiTCAD/verify/reference/special_functions.hpp`,
`src/verify/reference/special_functions.cpp` and `tests/verify/reference/special_functions_test.cpp`.

Every 15.11 check passes on its first run (Release):

| Check | Measured worst | Bound |
|---|---|---|
| erf(1), tanh(1) against OEIS | 0.50, 0.23 u_dd | 4 |
| Special values (li2 at −1, ½, 1; atan 1; asin ½; sin π/6; cos π/4) | 0.64 u_dd | 8 |
| erf against Gauss–Kronrod on (0, 6] | 2.5 u_dd | 16 + U_ref |
| erfc against Gauss–Kronrod on [3, 25] | 85 u_dd (within U_ref) | 16 + U_ref |
| atan on (0, 20]; li2 on [−20, 0.99] | 0.99; 1.5 u_dd | 16 + U_ref |
| log1p against Taylor | 0.43 u_dd | 4 |
| expm1∘log1p; log1p_exp antisymmetry; tanh doubling | 1.2; 0.11; 0.96 | 16; 8; 16 |
| sin² + cos²; sin∘asin; atan x + atan 1/x | each ≤ 1 u_dd | 8; 16; 8 |

Two notes on the build:
- **The erfc figure.** The 85 u_dd grows like x²/8. That is `exp`'s documented error at argument x² (TC1: "a few
  u_dd plus |x|/8 u_dd"), which both erfc and the integrand inherit. A diagnostic with the integral in 0.25-wide
  pieces agrees to 0.4 u_dd at x = 3 and 28 u_dd at x = 25. Nothing points at the continued fraction.
- **The domains of sin and cos** are |x| ≤ 2, not 15.11's ≤ 1, so that sin(asin x) can be formed up to π/2. No bound
  changed.

Mutations: 20 mutants of `special_functions.cpp`.
- 17 were caught on the first run. M16 (atan2's x < 0 branch), M17 (asin's (1 − x)(1 + x) written 1 − x²) and M20
  (erfc below 3) survived.
- Two checks were added (none loosened):
  - atan2 at three left-half-plane points, to 8 u_dd: 0.21 u_dd worst;
  - erfc on [0, 2.99] against Gauss–Kronrod, **absolutely**, to 16 u_dd + U_ref: 2.1 u_dd worst.
- After the additions, 19 of 20 are caught. M17 is equivalent for every input the references use: for a
  double-valued x the double-double x² is exact, so 1 − x² equals (1 − x)(1 + x).
- **Found by M20 and now documented in the header:** below x = 3, erfc is 1 − erf. It is accurate absolutely, not
  relatively; its relative error grows by erf/erfc, which is 4.5e4 at x = 3. No TVB-31 reference uses erfc there.
  MLDA's mean factor uses erf.

Tests: all 44 non-`full` verify ctests pass in Release, and `verify.special_functions` takes 0.18 s. The Debug cost
is measured at close-out (MI-Q9's 8 s allowance).

**15.13 Step 3 as built: TVB-31.01–31.05, the band set and narrowing (2026-11-06).**

Files:
- `verify/reference/model_band.hpp`: the source formulas, templated on the scalar (double-double for values, complex
  double for complex step). Varshni; the power-law N(T); Pässler as Couderc eq (10); Couderc's mass fits and N(T);
  ln n_i, n_i, E_c − E_i and d ln n_i/dT; Slotboom in both forms; n_ie. k is CODATA 2018's exact k/q.
- `verify/reference/decimal.hpp`: the register's printed decimals, read as the correctly rounded double and exactly
  as M·10^e in double-double.
- `verify/reference/condition.hpp`: S(q) = Σ|xᵢ ∂q/∂xᵢ| by double-double central differences, step 2⁻²⁷|xᵢ|.
- `verify/benchmark/tvb31_common.{hpp,cpp}`: 15.9's value rule, U_ref = 2⁻⁸⁰(|q_ref| + S) with TB-Q4, the register
  table check, and the register-built parameter sets of the formula check (15.2).
- `verify/benchmark/tvb31_band.cpp`: entries TVB-31.01–31.05.
- Tests: `tests/verify/reference/model_band_test.cpp` and `tests/verify/benchmark/tvb31_band_test.cpp`.
- The contract gains grade **L** (MI-Q5). It is treated as PC, with the reason "provenance insufficient (L: legacy,
  no source read)", and it is not a valid E grade. `contract_test` pins both.

L0 (`model_band_test`):
- Every analytic derivative the entries judge NiTCAD's partials against agrees with complex step to 64 ε S. These
  are dEg/dT for both gaps, d ln N/dT for both DOS forms, dM/dT, and d ln n_i/dT and d ln n_ie/dT for both sets. They
  are checked at 16 points a decade over 1–600 K. Each truncation estimate is below 1e-3 of its bound.
- Identities checked: Eg(0) = Eg0; Pässler's leading χ² term at 1 K; n_i² = Nc Nv e^(−Eg/kT); N(300 K) = N300; the
  documented Slotboom form equals the published one minus E0√½; k within CODATA's printed truncation.

Grids and judgements:
- **T grids:** 64 points a decade over the source's range (31.01 and the legacy half of 31.04: 1–500 K; 31.02,
  31.03 and the default half of 31.04: 1–400 K, Couderc's fits), plus 300 K. Two points beyond the range are
  reported `not_applicable`.
- **N grid:** 64 a decade over 1e14–1e21 cm⁻³, plus N0, its `nextafter` neighbours and N0(1 + 10⁻ʲ) for
  j = 1–12. The last set was added after mutant M-cut (below). Out of range: 1e13, 3e21, 1e22.
- **Cuts:** Slotboom's cut at N0 is judged by the value rule against the documented reference, which is 0 there, so
  NiTCAD must return exactly 0.
- **n_i** is judged where it is a normal double; below about 9.1 K it is `not_applicable`.
- **Couderc's printed n_i** (eq (15)) is a separate fit, so no printed value is judged in these entries.
- **Inputs:** the formula check runs NiTCAD's functions on the shipped sets with the register's printed values
  written over their fields. The reference reads those decimals exactly. The intrinsic level's depth takes the shipped
  χ (grade L) as an input.

Outcomes (Release, first run):

| Entry | Pass | Finding | Not applicable |
|---|---|---|---|
| 31.01 | every value and slope; 5 PB table rows | 2: eps_r and χ, grade L | 12 (beyond 500 K) |
| 31.02 | all; 4 PA rows | 0 | 4 |
| 31.03 | all; 8 PA rows | 0 | 16 |
| 31.04 | all but 2 | 2: n_i at 9.31 K, both sets, known finding (vi) | beyond range, and n_i below 9.1 K |
| 31.05 | the documented form throughout; 2 PB rows | the published form at every in-range point: known finding (F1) | 18 |

- **(vi), new (MI-Q12).** `physics::intrinsic_density` forms sqrt(Nc Nv)·e^(−Eg/2kT). Where Eg/2kT exceeds
  −ln DBL_MIN = 708.4 (about 9.1–9.6 K for silicon), the second factor is subnormal while n_i itself (about 1e-300)
  is still normal, so digits are lost. Measured: 1.2e3 and 6.8e3 of the bound at 9.31 K. Physically it does not
  matter, but it is a miss under the frozen rule.
  - **MI-Q12.** (A) Record (vi) as a known finding over the region where the factor is subnormal (`registry.cpp`
    computes the region from each set's own gap), and propose separately the fix exp(½ ln(Nc Nv) − Eg/2kT)
    (MI-Q10: no NiTCAD change in TVR-3). (B) Leave it a fail.
  - **Ruling (owner, "choose recommended", 2026-11-06): A.** The frozen criterion is unchanged; the test records the
    two findings and allows a pass at the region's edge (9.65 K passes).
- **(F1)** is confirmed as 15.10 recorded. The published form misses by E0√½ = 4.89 meV in ΔE_g, about 9.9 % in n_ie
  at 300 K, and E0√½/(2kT²) in the slope. The documented form passes at every point, including the cut.
- **The BGN values** match Klaassen et al. 1992 (PB), so their table rows pass. The header's attribution of them to
  Slotboom remains the F1 attribution item for the owner.

Mutations:
- **Reference** (`model_band.hpp`): 17 mutants.
  - 11 caught on the first run.
  - 5 survivors were below the criterion's 1e-13 floor (relative changes of 1e-15 to 1e-14). The same five at 1e-12
    (n_ie's exponent at 1e-10, where S is large) are all caught.
  - M6 (the Pässler exponent formed as the double 1/6) changes values by about 1e-17. It is a double-double guard the
    double-rounded criterion cannot see by design.
- **NiTCAD** (code under test, restored after each run): 16 mutants of `semiconductor.cpp`, `thermal.cpp`,
  `temperature.hpp` and `bandgap_narrowing.cpp`. These covered the gap, the DOS, n_i, the depth, the three slopes,
  the narrowing and its cut.
  - All are caught except M-cut, N0 moved to 1.001·N0. Between the cut's `nextafter` neighbour and the next grid
    point (3 % above N0) nothing was judged, so the points N0(1 + 10⁻ʲ) were added. With them, M-cut and a 1e-7 shift
    are both caught. Adding points tightens and widens nothing.
- **Baseline:** mutation baselines exclude nothing. The runs started only once (vi) was recorded and the suite was
  green; the first reference run, made against a red baseline, was discarded.

Cost: `[tvb31]` takes 0.24 s in Release. Debug is timed at close-out.

**15.14 Step 4 as built: TVB-31.06–31.10, mobility (2026-11-06).**

Files:
- `verify/reference/model_mobility.hpp`: Caughey–Thomas with the documented μ_max(T); Masetti (both sets);
  Canali and its right-hand slope at E = 0; Lombardi and its three partials; the energy-balance balance root (closed
  forms for β = 1 and 2, otherwise Newton), the owner's regularization, and the six returned quantities.
- `condition.hpp` gains a vector form, so one evaluation per perturbation serves several quantities.
- `verify/benchmark/tvb31_mobility.cpp`: entries TVB-31.06–31.10. The run functions take the tier, which sets grid
  density (`per_decade`, 64 in every tier until the close-out timing; MI-Q9).
- The register gains `laws`: the printed T-laws behind its law_300K rows. The table check judges a shipped value
  against the law at 300 K to the law's printed digits, i.e. the sum of the value's changes for half a unit in the
  last digit of the coefficient and of the exponent.
- Tests: `tests/verify/reference/model_mobility_test.cpp`, `tests/verify/benchmark/tvb31_mobility_test.cpp`, and the
  shared `tvb31_check.hpp`.

L0:
- Every analytic partial agrees with complex step: CT and Masetti dμ/dT, Canali dμ/dE, and Lombardi's
  dμ/dE, dμ/dT and dμ/dμ_b.
- 31.10's root is a double-double iteration, so its second method is Brent's zero (TC5). They agree to 1e-27
  relative over a = 1e-8–1e6 for β = 1, 1.5 and 2. The partials are checked against double-double central
  differences, and x and dx/da are continuous across a0.
- **MI-Q13 (an L0 reading, ruled "choose recommended" 2026-11-06, A).** 15.9 says each analytic derivative is checked
  "against complex step to 64 ε S". For dμ/dμ_b ≈ 1, S(d) is about 1e-4·|d|, so 64 ε S is below complex step's own
  rounding of about ε|d|. That check cannot be met by a correct derivative. (A) Apply 15.9's partial rule in full,
  max(1e-13 |d|, 64 ε S), which is "the same rule" as for values. (B) Leave the bare 64 ε S and report such points
  `not_evaluated`. **A.** It is applied in both L0 files (`model_band_test` had passed under the bare form).
  Nothing in the benchmark changes: there the rule was always the full one.

Outcomes (Release, first run):

| Entry | Pass | Finding | Not applicable |
|---|---|---|---|
| 31.06 | every value and slope; 2 PB T-exponent rows | 8: the grade-L CT rows (F3) | 12 |
| 31.07 | μ at 300 K and the documented slopes and T sweep; 17 PB rows | the published-form slope at every in-range N (F6); 2 tag-C T-exponent pins | 22 |
| 31.08 | everything with the shipped and the law parameters; the 2 v_sat rows | 2: β_n = 2 and β_p = 1 against the laws' 1.109 and 1.213 (F2) | 8 |
| 31.09 | all; 8 PB rows | 0 | 32 |
| 31.10 | all but (vii); 2 PA τ_w rows | the tag-C a0 pin; (vii) | 36 |

- **(F2), now judged.** The electron v_sat of 1.07e7 is within Jacoboni's law at 300 K (1.0705e7, ±3.4e5 for the
  law's printed digits). So is the hole v_sat 8.37e6, within ±2.7e5 of 8.344e6. The β values of 2 and 1 are outside
  Canali's β(T) at 300 K and are register findings. Both β sets pass the formula check.
- **(vii), new (MI-Q14).** `carrier_temperature_mobility` forms ∂μ/∂μ0 = m + m_a·a and ∂E_eff/∂μ0 =
  v_sat(x_a·a − x)/μ0². Both cancel as a grows: each tends to 0 while its terms stay O(m). From a ≈ 25 (β = 2) up to
  a ≈ 1100 (β = 1), and through the range above that, the error is up to 1e8 times the bound. It is about ε·m in
  absolute terms, so harmless beside the Jacobian's O(m) entries, but a miss of a well-conditioned partial.
  - **MI-Q14.** (A) Record (vii) as a known finding for a ≥ 20 on those two partials, and propose separately the
    cancellation-free forms (MI-Q10). (B) Leave them failing.
  - **Ruling (owner, "choose recommended", 2026-11-06): A.** The test allows a pass inside the region.
- **Register corrections, made before any outcome depended on them.**
  - 31.08's electron law value was typed 1.0704e7 in step 1. The arithmetic gives 1.53e9·300^−0.87 = 1.07052e7, so
    the row is now 1.0705e7. The table check refuses (`not_evaluated`) a law_300K row that does not round its law,
    which is how it was found. The F2 bullet of 15.10 should read 1.0705e7.
  - 31.10's tag-C row named `energy_balance.canali_beta`, which binds no NiTCAD field. The ruling it cites (Unit 52,
    "β = 1") fixes the regularization threshold a0 = 1e-8, and the row is rebound to that.

Mutations:
- 13 of `model_mobility.hpp`: all caught.
- 13 of NiTCAD's `mobility.cpp`, `thermal.cpp` (the CT slope), `field_mobility.cpp`, `surface_mobility.cpp` and
  `energy_balance.cpp`. 12 caught on the first run. The survivor moved the a ≤ 0 cut to a ≤ 1e-30: the positive grid
  started at 1e-12. Points a = 1e-20, 1e-40 and 1e-100 were added, and it and a 1e-90 variant are caught.

Cost: 31.10 took 11.8 s in Release on the first run, from re-solving the root once per quantity and per
perturbation. With the vector condition sums and a double Newton start before the double-double steps it takes
0.61 s. All of `[tvb31]` takes 1.65 s in Release. That puts the Debug figure above the 8 s allowance (MI-Q9); the
fast/full split is settled at close-out against the measured Debug time.

**15.15 Step 5 as built: TVB-31.11–31.16, recombination and traps (2026-11-06).**

Files:
- `verify/reference/model_recombination.hpp`:
  - the Scharfetter lifetime;
  - SRH, Auger and radiative, each with dR/dn and dR/dp and the chain through E (E's partials are inputs);
  - the trap kinetics of one level (f, 1 − f, U/N_t and four partials);
  - the equilibrium occupancy;
  - the trap-band density in closed form, D kT[ln(1 + e^−u_lo) − ln(1 + e^−u_hi)], and by adaptive Gauss–Kronrod (TC3).
- `verify/benchmark/tvb31_recombination.cpp`: entries TVB-31.11–31.16. `bound_field` gains the lifetime and Auger
  fields.
- Tests: `tests/verify/reference/model_recombination_test.cpp`, `tests/verify/benchmark/tvb31_recombination_test.cpp`.

Grids:
- 31.12–31.15: n on the 1e-10–1e21 log grid (1985 points), with p = r·E/n for r ∈ {1e-6, 1, 1e6}. r = 1 gives the
  zero. Points with p outside the range are `not_applicable`.
- 31.12 and 31.13: the Boltzmann E and the Fermi–Dirac E of `statistics.hpp`.
- 31.14: the four shipped non-zero B values. GaAs is also run with the Fermi–Dirac E, added after mutation (below).
- 31.15:
  - `fermi_occupancy` over x ∈ [−750, 750] in steps of ½, plus 0 and its neighbours; |x| > 745 is
    `not_applicable`;
  - the kinetics at levels τ ∈ {−20, 0, 20} × r, plus one level with dn1/dn and dp1/dp non-zero.
- 31.16: four bands up to the gap, both densities. E_F runs −0.5..0.5 eV in 1/256 at 300 K, and T runs 200–500 K at
  E_F 0 and 0.25 eV.

L0:
- **A defect in the first reference, fixed before any benchmark run.** The partials were first written by the bare
  quotient rule, (A·D − B·C)/D², and that form cancels when one capture term dominates D. At n = 1e-10, p = 1e36 the
  SRH terms agree to 1e-26 relative, below even double-double. The reference now uses the expanded forms, in which the
  cancellation is symbolic:
  - SRH: dR/dn = (τ_n p(p + n_ie) + τ_p(p n_ie + E) − E_n D)/D²;
  - trap: df/dn = c_n(emp − n1′·occ)/D², and dU/dn = c_p(p + p1)·df/dn;
  - and the mirrors of both for p.
- Complex step is forward-mode differentiation of the formula as written, so it inherits the same cancellation. The
  second method therefore differentiates algebraically equal value forms that are free of it:
  - SRH as n/τ_n − [n τ_p(n + n_ie)/τ_n + n n_ie + E]/D, or its mirror, chosen by real parts;
  - the smaller of f and 1 − f;
  - U in its hole or electron capture form.
- h is scaled to the input (1e-30|x|). At the default 1e-200 these small rates underflow in the imaginary part.
- SRH agrees with the Brent root (TC5) of the capture–emission balance. The trap U agrees with c_n(n(1 − f) − n1·f).
  f + (1 − f) = 1. The band closed form agrees with Gauss–Kronrod to 1e-26 relative plus the quadrature's estimate.

Outcomes (Release):

| Entry | Pass | Finding | Not applicable |
|---|---|---|---|
| 31.11 | 1028 | 3: τ_n0 and τ_p0 (L), N_ref (PC) | 2 |
| 31.12 | 23646 | 1644: (viii) | 10440 |
| 31.13 | 25292, incl. the 2 PB Auger rows (bit-equal) | 0 | 10440 |
| 31.14 | 63225 | 0 | 26100 |
| 31.15 | 101277 | 6305: (ix) | 40380 |
| 31.16 | 2488: the 6-point Gauss–Legendre panels ≤ kT meet the value rule everywhere | 0 | 48 |

- **(viii) and (ix), new (MI-Q15).** NiTCAD's `srh_recombination` and `trap_kinetics` form their partials by the bare
  quotient rule, the form the reference had to drop. The miss is up to about 800 times the bound:
  - SRH, at np/E = 1e6 from n ≈ 1.2e5, and marginally (1.3–1.7 times the bound) with the degenerate Fermi–Dirac E at
    n ≥ 9e20;
  - the trap partials, wherever f or 1 − f is near 1 or one capture term dominates.
  The absolute errors are rounding of the shipped form, about ε times its terms.
  - **MI-Q15.** (A) Record both as known findings, and propose the expanded forms separately (MI-Q10). (B) Leave
    them failing.
  - **Ruling (owner standing instruction, "choose recommended", 2026-11-06): A.**
  - **Region, as built.** The first version used parameter regions: n ≥ 1e5 for SRH, and everywhere for the trap
    partials. Mutation showed that "everywhere" absorbed real defects: a dropped dp1/dp factor became a finding. The
    finding is now defined by its mechanism, inside the entry:
    - a partial that misses the value rule is (viii) or (ix) only where |q_sim − q_ref| ≤ 64 ε T, with T =
      (|A|·D + |B·C|)/D² the size of the shipped form's terms;
    - a larger miss, such as a wrong term (of order T), stays a fail.
    Every miss of the first run is within the floor. The pass criterion is unchanged (TVM2). The registry's
    known-finding table stays at 25 entries, and the entries emit a stale note if their finding stops reproducing.

Mutations:
- 15 of `model_recombination.hpp`. 12 were caught on the first run. Two survivors were closed by new L0 checks: the
  Scharfetter lifetime (exact at N = 0 and N_ref), and the trap U against its capture form. A sub-floor ×(1 + 1e-20)
  on the band closed form is caught at 1e-12.
- 19 of NiTCAD's `recombination.hpp/.cpp` and `interface_traps.hpp/.cpp`. 14 were caught on the first run (one by a
  build failure, rerun in a compilable form and caught). The survivors:
  - B(p − E_n) → B·p: 31.14 had only the Boltzmann E. GaAs's Fermi–Dirac case was added, and it is caught.
  - Three mutants of the dn1/dn and dp1/dp terms were absorbed by the "everywhere" region. They are caught under the
    floor rule above.
  - A 1e-12 shift of the Gauss–Legendre nodes was sub-floor. It is caught at 1e-9.

Cost (Release): 31.11–31.14 take 0.5 s, 31.15 0.4 s and 31.16 1.1 s, about 2 s in all. The fast/full split
(MI-Q9) is settled at close-out.

**15.16 Step 6 as built: TVB-31.17–31.19, ionization and the neutral equilibria (2026-11-06).**

Files:
- `verify/reference/model_statistics.hpp`:
  - the ionized density, with the bound fraction and the documented complete ionization at depth ≤ 0;
  - Altermatt–Schenk E_dop(N) and b(N);
  - the neutral equilibrium as the Brent root (TC5) in η of the double-double charge balance. F_1/2 comes from the
    TVB-02 reference: its series for x ≤ −2, tanh–sinh (TC4) above.
  - The condition sums of n, p and η at the root are by implicit differentiation, dη/dx_i = −(∂G/∂x_i)/(∂G/∂η).
    The densities are held fixed for the inputs that do not move them; they are linear in n_ie, and only g_n and g_p
    are differenced through F_1/2.
- `verify/benchmark/tvb31_statistics.cpp`: entries TVB-31.17–31.19. `bound_field` gains the ionization fields, and
  `other_shipped_value` the Altermatt sets. The rounding-floor findings of 15.15 (`FloorFinding`, `classify_floor`)
  move to `tvb31_common.hpp`.
- Tests: `tests/verify/reference/model_statistics_test.cpp`, `tests/verify/benchmark/tvb31_statistics_test.cpp`.

Grids:
- 31.17 (N = 1e17, kT at 300 K): η −40..40 in steps of ¼ at depths 0, 0+ (the next double), 0.045, 0.1 and 0.3 eV;
  depth 0..0.3 eV in 256 steps at η −40, −5, 0, 5 and 40; g 2 and 4; b 1, 0.5 and 1e-3.
- 31.18: N 0 and 1e14–1e21, with the register's Table I values written into the function's argument.
- 31.19 (legacy silicon at 300 K):
  - complete ionization, with both statistics, at net doping 0 and ±1–1e21;
  - ionized, with both statistics and both the legacy fixed levels and the Altermatt levels, over N_D or N_A
    1e10–1e21, plus 25 compensated pairs N_D = N_A(1 ± 10^-k) about N_A = 1e17.

L0:
- dN_ion/dη agrees with complex step (h = 1e-30|η|). At b = 1 the value is N·s exactly, and complete ionization holds
  at depth 0.
- Altermatt is exact at N_ref and N_b (c = 2).
- The F_1/2 series and Gauss–Kronrod agree at the switch, and tanh–sinh agrees with Gauss–Kronrod from −1.999 to 60,
  to 1e-26.
- The Boltzmann root agrees with its closed form ln((C + √(C² + 4n_ie²))/(2n_ie)) to 1e-29. The Fermi–Dirac root
  satisfies n − p = C, and degeneracy raises |η| at |C| ≥ 1e19.
- The implicit condition sums agree with `condition_sum` over re-solved roots to 1e-6.

Outcomes (Release):

| Entry | Pass | Finding | Not applicable |
|---|---|---|---|
| 31.17 | 34680 | 4: the legacy levels and degeneracies (PC) | 180 |
| 31.18 | 1810, incl. the 10 PA Table I rows (bit-equal) | 0 | 4 |
| 31.19 | 32424: every n and p, and every η but (x) | 942: (x) | 12 |

- **(x), new (MI-Q16).** The η of `fermi_dirac_neutral_equilibrium` misses where |η| is small (|C| ≲ 6e7, η ≲ 3e-3).
  - Mechanism: the residual forms ln c + log1p(minor/c), with c = |C|/n_ie. Its two terms are each about |ln c| and
    cancel when c ≪ 1, which leaves an absolute error of about 5e-15. The header documents the root as "good to a
    few 1e-14 in y".
  - n and p pass, because the majority is formed as |C| + minor.
  - **MI-Q16.** (A) Record (x) as a known finding. It is classified by the residual's rounding floor, T = |η| +
    2|ln c|, as for (viii) and (ix), and a fix is proposed separately (MI-Q10). (B) Leave it failing.
  - **Ruling (owner standing instruction, "choose recommended", 2026-11-06): A.**
- **31.19's first run took 151 s in Release.** It used Gauss–Kronrod F_1/2 (1.4 ms each), a ±60 Brent bracket, and
  condition sums that re-evaluated both densities for all 11 inputs. With tanh–sinh above x = −2 (0.3–1.1 ms), the
  bracket started about NiTCAD's η (its ends verified by sign and widened ×4, so the root is still the reference's
  own balance), and the fixed-density condition sums, it takes 6.3 s. That is too long for `fast` at 64 points a
  decade, and is for the close-out split (MI-Q9).

Mutations:
- 8 of NiTCAD's `ionization.cpp`, 4 of `statistics.cpp` (the Fermi–Dirac neutral solve), and 8 of
  `model_statistics.hpp`.
- **The first run was invalid.** Its Catch filter, `TVB-31.1[789]*`, matched no test, because Catch has no
  character classes. "No tests ran" exits non-zero, so every mutant looked caught. Rerun with
  `TVB-31.17*,TVB-31.18*,TVB-31.19*`, after a green baseline check: 19 of 20 are caught. The survivor writes s =
  1 − e·d for d on y ≤ 0, which is the same number to rounding (an equivalent mutant). Step 5's runs used
  `TVB-31.1*`, which matches the ten entries 31.10–31.19, and stand.

**15.17 Step 7 as built: TVB-31.20–31.22, impact ionization and band-to-band generation (2026-11-06).**

Files:
- `verify/reference/model_generation.hpp`:
  - vOdM α = γA exp(−γB/E) with dα/dE; γ(T) = tanh(h/300)/tanh(h/T), h = ħω/2k, with dγ/dT; Kane G = AF² exp(−B/F)
    with dG/dF.
  - The direct-gap WKB functions of Esseni et al. (2017):
    - the edge-safe α and cos θ (1 ± α by their rationalized forms);
    - κ (eq 9), and its antiderivatives in δ for κ and 1/κ;
    - the segment integrals with Leibniz partials.
  - Segment brackets:
    - For D ≥ 1e-2 the brackets D f(v_max) − ∫f are closed form.
    - Below 1e-2 they are tanh–sinh of the difference.
    - Where the in-gap part lies within 1e-3 of its distance to the nearer edge, one Kronrod-15 rule computes both the
      brackets and the integrals. That interior rule was added after the first run, when tanh–sinh's estimate stalled
      at D = 1e-13 (a reference defect, not a NiTCAD one); the L0 checks it against tanh–sinh at D = 1e-5 and against
      the D → 0 limit.
  - The path rate (eq 11) and its partials; dG/dI_ik takes the bracket 1 − e^(−z)(1 + z) by its series for z < 1.
  - eq (8).
- `verify/benchmark/tvb31_generation.cpp`: entries TVB-31.20–31.22.
  - `bound_field` gains the vOdM coefficients, the phonon energy and the Kane pair.
  - `other_shipped_value` gains `impact_ionization.gamma_at_300K`, NiTCAD's γ at 300 K.
  - `tvb31_common.hpp` gains `cut_record`, 15.9's documented-cut rule: q_sim = 0 and |q_ref| ≤ |q_ref at the
    threshold|.
- The register's 31.20 tag-C row is rebound to `impact_ionization.gamma_at_300K = 1`. It pins T₀ = 300 K through its
  consequence: γ(300 K) = tanh(x)/tanh(x) is exactly 1. The row replaces the unbindable "T₀" of 15.10; the form stays
  F4's.
- Tests: `tests/verify/reference/model_generation_test.cpp` (412 assertions) and
  `tests/verify/benchmark/tvb31_generation_test.cpp`. The registry holds 32 benchmarks.

Grids:
- 31.20, legacy silicon with the register's printed coefficients:
  - γ and dγ/dT at T 200–500 K on 64 points a decade, plus 300 K, and 150 and 600 K beyond.
  - α and dα/dE for both carriers at T 200, 300 and 500 K, over:
    - E 0 and 1e4–1e7 V/cm, and 1e8 beyond;
    - the switch field and its neighbouring doubles;
    - the cut field γB/700 with its neighbours and γB/700·(1 ± 10^-j), j = 1–12.
  - The cut and switch points are judged wherever they fall.
- 31.21: the same field grid with the shipped Kane pair (grade L) and the cut B/700.
- 31.22, two test bands (no shipped material has direct-gap masses): u = 14.9 (E_g 1.42 eV, m_r 0.0335 m₀) and u = 1.5
  (0.5 eV, m₀/3).
  - κ and both antiderivatives at δ in steps of 1/1024, at 10^-k and 1 − 10^-k (k = 1–15), and at −0.5 and 1.5.
  - Segments with L = 1e-8 m: 90 ordered pairs of ten endpoints in and around the gap, and 96 near-coincident pairs
    v, v + 10^-k (k = 3–14) about v = 1e-6, 0.3, 0.6 and 0.999999, in both directions.
  - eq (8) at F 1e7–1e10 V/m, and at 1e6, where the exponent exceeds 700.
- The path rate at 672 points: F ±1e7–1e9 V/m, I_k 0–340, I_ik 0 and 1e-20–1e-14 m², and z = k_m² I_ik from 1e-10 to
  800.

Outcomes (Release, 0.6 s for the three):

| Entry | Pass | Finding | Not applicable |
|---|---|---|---|
| 31.20 | 2745, incl. the PB phonon-energy row and every cut point | 8: the seven vOdM coefficients (PC, F4) and the tag-C γ(300 K) pin (I0) | 16 |
| 31.21 | 442, incl. every cut point | 2: the Kane pair (grade L, F5) | 2 |
| 31.22 | 11314 | 994: (xi) 866, (xii) 128 | 14 |

- **MI-Q17 (an entry grade).** 31.20's form was first written as PC, which would have made every α comparison a
  finding. The form α = A exp(−B/E) is the one printed by Allpix Squared and by Sentaurus's summary (two
  reproductions, so PB); only the coefficients are PC (one reproduction of the abstract, F4). γ(T) is a NiTCAD choice
  (tag C), so its comparisons are PD: its documented form, checked as written.
  - (A) PB for the form, PD for γ's comparisons, with the pins on the register rows. (B) PC for the whole entry.
  - **Ruling (owner standing instruction, "choose recommended", 2026-11-06): A.** It was ruled before the first run.
- **(xi), new (MI-Q18).** `wkb_segment` forms I = L(K(hi) − K(lo))/D and its partials L(f(v) ∓ dK/D)/D.
  - Mechanism: K is about u·π/4 times scale/(2uħ) at each end. A short segment therefore loses ε·K/(κD) relative in
    the value and ε·K/(κD²) in the partials. Example: I_k over [0, 1e-9] is off by 0.3% (9.8e-7 absolute, against
    64ε·T = 1.4e-6).
  - Below D = 1e-10 the documented midpoint rule takes over. It truncates L f″D²/24 in the value and drops the
    ±L f″D/12 of each partial, which is 2e-6 relative at δ = 1e-6, D = 1e-11.
  - Classification:
    - for D ≥ 1e-10, the rounding floor with T = L(K_t(hi) + K_t(lo))/D (the partials add L·f(v)/D and divide by D),
      where K_t and G_t are the sizes of the antiderivatives' terms;
    - below 1e-10, twice the midpoint rule's leading truncation plus 64ε|q_ref|, with f″ by double-double central
      differences.
  - Every miss falls within one of the two bounds.
- **(xii), new (MI-Q18).** `wkb_path_rate`'s dG/dI_ik = pre·t·(k_m² e/I_ik − (1 − e)/I_ik²) cancels as z → 0
  (relative ε/z²). Floor: T = pre·t·(k_m² e (1 + z)/I_ik + (1 − e)/I_ik²). All 128 misses (z = 1e-10 and 1e-4) are
  within it.
- Fixes proposed separately (MI-Q10):
  - for (xi), the brackets by quadrature or Taylor series of the difference below D ≈ 1e-2, and the midpoint partials
    with their ±f″D/12 term;
  - for (xii), the series of 1 − e^(−z)(1 + z) for z < 1, as the reference uses.
  - The practical effect is small: mesh segments in δ are rarely below 1e-4, and z ≪ 1 means a path with almost no
    final states.
- **MI-Q18.** (A) Record (xi) and (xii) as known findings, classified as above, and propose the fixes separately. (B)
  Leave them failing.
  - **Ruling (owner standing instruction, "choose recommended", 2026-11-06): A.**
  - The frozen 31.22 test, written before the first run, expected no findings. It now names the two regions; the
    criterion is unchanged (TVM2).
- **31.20 and 31.21 met their frozen outcomes on the first run.**

Mutations, each run after a green baseline with the exact filter `TVB-31.20*,TVB-31.21*,TVB-31.22*`; all 23 are
caught:
- 15 of NiTCAD's code:
  - `impact_ionization.cpp`: the switch's ≥ → >, the cut at 699, γ dropped from α, dα/dE · (1 + 1e-12), and
    T₀ = 300.0000001;
  - `thermal.cpp`: the slope zeroed above 400 K;
  - `band_to_band.cpp`:
    - Kane: dG/dF's x · (1 + 1e-12) and the cut at 701;
    - the WKB functions: S's ¼ + 1e-10, K's c³/3.0000001, G's 0.9999999·c, the midpoint switch at D < 1e-4, in_lo
      with ≥ 0, I_k · (1 + 1e-12) in the path rate, and eq (8)'s E_g^1.5000000001.
- 8 of `model_generation.hpp`, caught by its L0: dα/dE, dγ/dT, sinh → cosh, both antiderivatives, the g(z) series,
  eq (8)'s 18π², and the Kronrod bracket · (1 + 1e-9).

**15.18 Step 8 as built: TVB-31.23–31.29, tunnelling, emission, MLDA and thermal (2026-11-06).**

Files:
- `verify/reference/model_tunnelling.hpp`:
  - Gehring's segment action L·g(a, b) with its partials (rationalized where both ends are inside, the
    turning-point form otherwise);
  - the exponent S(E), dS/dE and TC = e^(−2S);
  - the supply N(E) with its logistic partials;
  - the current, as an adaptive Gauss–Kronrod integral between every barrier end inside NiTCAD's cut range.
- `verify/reference/model_thermal.hpp`:
  - the emission velocity in both forms, with dv/dT through N(T);
  - MLDA's λ, f and the mean factor (its series below a/λ = ¼);
  - the conductivity power law;
  - F₋₁, F₀ and F₁ (Li₂'s inversion above 0; e^x scaled by 2^512 below x = −600, outside double-double's domain);
  - h, and dh/dx in a cancellation-free form, 2(S₀² − S₁S₋₁)/S₀². The numerator is summed as Σcₙ(−z)ⁿ with
    constant cₙ = Σᵢ[1/((i+1)(n−i+1)) − 1/(i+1)²], so dh/dx → e^x/2 without cancellation. The first draft took
    2 − 2·product in double-double, which is itself only rounding below x ≈ −69; this was found by analysis
    before the first run.
- `verify/benchmark/tvb31_tunnelling.cpp` (31.23–31.25) and `tvb31_thermal.cpp` (31.26–31.29).
  - `bound_field` gains the three thermal fields, and `other_shipped_value` the two MLDA masses.
  - The check helper `tvb31_check.hpp` accepts a named set of `not_evaluated` comparisons (TB-Q4).
- Tests: `model_tunnelling_test.cpp` and `model_thermal_test.cpp` (207 assertions together, after the mutation runs),
  `tvb31_tunnelling_test.cpp`, `tvb31_thermal_test.cpp`. The registry holds 39 benchmarks.
- **G-tag correction.** Steps 5–7 copied mobility's G6 into every entry. They now carry their own governing
  equation: recombination G7, statistics G4, generation G8 (31.22 G8 and G12). This is a label only; no outcome
  moves.

Grids:
- 31.23:
  - L·g and both partials at the 13 × 13 pairs of a, b ∈ {−1, −1e-3, −1e-12, 0, 1e-12, 1e-6, 1e-3, 0.1, 0.5, 1, 2,
    3.1, 5} eV, and at a(1 ± 10^-k) about a = 1e-3, 0.5 and 3.1 eV (k = 1–15), with L = 1.5 nm;
  - √(2m₀q)/ħ;
  - S, dS/dE and TC of four barriers over E −1 to 3.5 eV in 1/64 eV, plus every barrier end and its neighbours:
    rectangular 3.1 eV, trapezoidal 3.1 → 2.1 eV, Fowler–Nordheim triangular 3.1 → −0.9 eV (all 1.5 nm, 0.5 m₀),
    and a two-layer 1 nm (0.5 m₀) + 2 nm (0.18 m₀) stack.
- 31.24: N and its three partials at 300 K, with x₁ = (F₁ − E)/kT from −64 to 64 in steps of ⅛, for F₁ − F₂ ∈ {0,
  1e-9, 1e-3, 0.1, 1} eV. |x| > 60 for either level is not applicable.
- 31.25: the four barriers at 300 and 400 K, for F₁ − F₂ ∈ {0.01, 0.1, 0.5, 1} eV with E_low = 0.
  - None of those barriers has an end inside the cut range. After the first run, one barrier whose ends (0.8 and
    0.3 eV) lie inside it was added, to exercise MI-Q8's kinks. These are diagnostics only, so no criterion moves.
- 31.26:
  - v(N, T) at N 1e15–1e21 and T 200, 300 and 500 K;
  - v and dv/dT of both carriers at T 200–500 K (150 and 600 K beyond), on the legacy and Couderc sets, each with
    the density-of-states form and with test Richardson constants of 120 A/(cm² K²) (no shipped set has any);
  - the Schottky overload's v.
- 31.27:
  - λ for both masses at T 200–500 K;
  - f at z/λ 0–12 in steps of 1/64;
  - the mean factor at a/λ 1e-8–1e3, with the switch ½ and its neighbours.
- 31.28: κ and dκ/dT on the register's printed values (1.48 W/(cm K), 1.3) at T 200–600 K, and the three register
  rows.
- 31.29: F₋₁, F₀, F₁, h and dh/dx at x −750..750 in steps of ¼, at −2..2 in 1/1024, and at the switch points ±1 with
  their neighbours; ±760 lie beyond the range.

Outcomes (Release, 0.9 s for the seven):

| Entry | Pass | Finding | Not evaluated | Not applicable |
|---|---|---|---|---|
| 31.23 | 4471 (every comparison) | 0 | 0 | 0 |
| 31.24 | 16553 | 1295: (xiv) | 0 | 2652 |
| 31.25 | — (41 diagnostic notes, MI-Q8) | — | — | — |
| 31.26 | 1603 | 0 | 0 | 92 |
| 31.27 | 1405 | 2: the PC quantization masses | 0 | 130 |
| 31.28 | 66 | 3: 1.48 (PC), the exponent's mismatch (1.33 against the printed 1.3), the heat capacity (L) | 0 | 4 |
| 31.29 | 47166 | 2828: (xiii) | 516: the underflow tail (MI-Q19) | 10 |

- **31.25 diagnostics (MI-Q8, for the owner's threshold).** NiTCAD's composite 4-point Gauss–Legendre agrees with
  the double-double integral:
  - to 1.1e-13 relative or better, wherever no barrier end lies inside the range (32 points);
  - to 1.2e-9 at 300 K and 2.8e-8 at 400 K, on the kinked barrier (the turning point crossing a segment end).
- **MI-Q19 (TB-Q4 below 2^-969, ruled before the first run).**
  - The rule: where |q_ref| < 2^-969, the double-double's low part is subnormal, so its uncertainty is at least the
    subnormal spacing 2^-1074. It is not 2^-80|q|, and a subnormal q_ref's own rounding is of the same size.
  - (A) Floor U_ref at 2^-1074 there, in the shared `value_record`. TB-Q4 then makes such a point not_evaluated
    wherever that exceeds 1e-3 of the bound. (B) Leave subnormal points judged by the value rule, which no double
    can meet below the subnormal spacing.
  - **Ruling (owner standing instruction, "choose recommended", 2026-11-06): A.**
  - Effect beyond step 8: 260 points of 31.15's occupancy tail at η = −745 (f_eq and df_eq/dη, about 1e-313 to
    5e-324) move from pass (bit-equal) to not_evaluated. Its frozen test now names them, and its other outcomes are
    unchanged (101017 pass, 6305 (ix), 40380 not applicable). No other entry moves.
  - By TB-Q4, an entry with any not_evaluated point is not_evaluated overall: 31.15 and 31.29 report so. `nitcad
    verify` does not fail on it.
- **MI-Q20 (an entry grade).** 31.28's power law is the legacy's documented form. One reproduction (Palankovski) is
  PC. (A) Judge the formula as PD (the documented form, checked as written, as MI-Q17's γ), with the grades on the
  register rows. (B) PC for the whole entry.
  - **Ruling (owner standing instruction, "choose recommended", 2026-11-06): A.** It was ruled before the first run.
- **(xiii), new (MI-Q21, ruled before the first run, from the shipped form).** `fermi_dirac_thermal_diffusion`'s
  dh/dx = 2 − 2F₁F₋₁/F₀² cancels as x → −∞, where dh/dx → e^x/2. The absolute error is about ε, so the relative
  error exceeds 1e-13 below x ≈ −6.
  - Classification: the rounding floor T = 4 (two terms of size 2).
  - All 2828 misses lie within it, between x ≈ −6 and the underflow tail.
  - Fix, proposed separately (MI-Q10): the series form above, or 2(F₀² − F₁F₋₁)/F₀² with its leading cancellation
    removed. The effect is an absolute 1e-16 in a thermal diffusion factor of ½.
  - **Ruling: A** (record and classify), as for (viii)–(xii).
- **(xiv), new (MI-Q22, ruled after the first run found it).** `tsu_esaki_supply`'s dN/dE = σ(x₂) − σ(x₁) cancels
  where both levels are many kT above E (both σ → 1).
  - The misses are at x₁ ≥ 7.75. For example, at F₁ − F₂ = 1e-9 the shipped value is exactly 0 against −2.9e-17.
  - Classification: the floor T = σ(x₁) + σ(x₂). All 1295 misses lie within it.
  - Fix, proposed separately: (e^(−x₁) − e^(−x₂))/((1 + e^(−x₁))(1 + e^(−x₂))) through expm1. The effect is an
    absolute ε in the current's energy partials.
  - **Ruling (owner standing instruction, "choose recommended", 2026-11-06): A.** The frozen 31.24 test now names
    the region; the criterion is unchanged (TVM2).

Mutations (after a green baseline; filter `TVB-31.23*,...,TVB-31.29*` on `nitcad_verify_benchmark_test`). There
were 18 one-literal NiTCAD mutants, each reverted by `git checkout`, and **all 18 were caught**:
- T01–T08, `gate_tunnelling.cpp`:
  - both inside partials;
  - the outside g;
  - the turning-point partials;
  - the scale (1e-10);
  - dS/dE;
  - the supply's dN/dF₂;
  - the log1p;
- T09–T10, `thermionic_emission.cpp`: the density-of-states form; the Richardson form;
- T11, T12, T16–T18, `thermal.cpp`:
  - dv/dT's ⅓;
  - the fit-less ½;
  - dκ/dT;
  - F₁'s π²/6;
  - h's 1.5;
- T13–T15, `quantum_confinement.cpp`: λ's 100; the mean factor's switch ½ → 0.01; its series truncated at 6
  terms.

`git status src/physics` is clean afterwards.

Reference mutants, against the L0 `[model_tunnelling],[model_thermal]`: 10, restored from a backup copy.
- In the first run, 8 were caught. Two survived: R01 (the scale's 100 → 100.0000001) and R03 (TC's −2 →
  −2.0000001). Nothing in the L0 pinned either.
- An L0 case was added: the scale against another arrangement of its constants, and a rectangular barrier's S and
  TC against κL√(V − E) and e^(−2S). Both mutants are now caught.
- The L0 is 207 assertions in 8 cases.

**15.19 Step 9: material tables, sources, criteria and rulings (frozen 2026-11-06, before code).**

Sources read for 31.31 and the radiative coefficients:
- The Ioffe NSM archive pages (http://www.ioffe.ru/SVA/NSM/Semicond/<material>/), read 2026-11-06 over plain
  HTTP; the fetch tool refused the host's certificate chain.
  - Ge, GaAs and GaInAs: `basic`, `bandstr`, `electric` and `thermal`.
  - SiC: the same four, plus `ebasic` and `recombination`.
  - AlGaAs: `basic`, `bandstr`, `ebasic` and `recomb`.
- The values are typed into the register as printed. Highlights:

  | Material | Printed values |
  |---|---|
  | Ge | ε 16.2; E_g = 0.742 − 4.8e-4 T²/(T + 235); N_c 1.0e19; N_v 5.0e18; χ 4.0; μ ≤ 3900 and ≤ 1900; μ_n ≈ 4.9e7 T^-1.66 and μ_p ≈ 1.05e9 T^-2.33; B 6.41e-14; Auger ~1e-30; P 0.013 and B 0.011 eV |
  | GaAs | ε 12.9; E_g = 1.519 − 5.405e-4 T²/(T + 204); N_c 4.7e17; N_v 9.0e18; χ 4.07; μ ≤ 8500 and ≤ 400; B 7.2e-10 at 300 K (Varshni 1967); Auger ~1e-30; Si donor ~0.006 and C acceptor ~0.02 eV (Milnes 1973) |
  | Ga₀.₄₇In₀.₅₃As | ε 13.9; E_g 0.74 eV (300 K), with Paul et al.'s E_g(x, T) and no Varshni pair; N_c 2.1e17; N_v 7.7e18; χ 4.5; μ < 12·10³ and < 300; B 0.96e-10; Auger 7e-29 |
  | 4H-SiC | No ε of its own ("the value of 6H-SiC is usually used"); E_g 3.23 eV (300 K) and E_g(0) − 6.5e-4 T²/(T + 1300); N_c ≈ 3.25e15 T^(3/2) and N_v ≈ 4.8e15 T^(3/2); μ ≤ 900 and ≤ 120; longest lifetimes τ_p 6.0e-7 s (n-type) and τ_n 1e-9 s (p-type); B 1.5e-12; C_n 5.0e-31 and C_p 2.0e-31 (Galeskas 1997); κ 3.7 W/(cm K); density 3.211 g/cm³; specific heat 0.69 J/(g K) (6H, "usually used") |
  | Al_xGa₁₋ₓAs | ε 12.90 − 2.84x; E_g 1.424 + 1.247x (x < 0.45, 300 K); N_c 2.5e19 (0.063 + 0.083x)^(3/2) (x < 0.41); N_v 2.5e19 (0.51 + 0.25x)^(3/2); χ 4.07 − 1.1x (x < 0.45); μ_n 8·10³ − 2.2·10⁴x + 10⁴x² (0 < x < 0.45); μ_p 370 − 970x + 740x²; B ~1.8e-10 |

- Not printed there: Ge and GaAs saturation velocities (figures only); SiC's χ and its N and Al levels (the code
  cites Ayalew's thesis, not read); any Canali β. These rows are grade L.

Questions (owner standing instruction, "choose recommended", 2026-11-06; each ruled A before the first run):
- **MI-Q23 (composition laws).** AlGaAs's fields are functions of x, and Ioffe prints them as laws in x.
  - (A) Judge NiTCAD's `algaas_parameters(x)` at x = 0, 0.05, …, 0.45 against each printed law, evaluated in
    double-double. The tolerance is the law's change for half a unit in the last printed digit of each coefficient,
    summed: law_check's NPARC rule (11.5) carried to a polynomial, and for the (a + bx)^(3/2) forms by their
    first-order sensitivity. Points beyond a law's printed range are not_applicable. The 300 K gap is judged as
    NiTCAD's E_g(300 K) of the set.
  - (B) Register the coefficients and compare them bit-equal. NiTCAD stores no coefficients, only the function.
  - **Ruling: A.**
- **MI-Q24 (grades of the Ioffe values).**
  - (A) PA where NiTCAD's code cites the Ioffe database for that value, as it does for the radiative coefficients of
    Ge, InGaAs, 4H-SiC and AlGaAs. Otherwise PC: one compilation, against the code's "handbook (Adachi)" or another
    unread citation.
  - A value printed with a qualifier is typed as its printed number, with the qualifier in the row's note: ~, ≤ and
    < (a mobility bound is the lattice mobility, μ_max).
  - (B) PC for all.
  - **Ruling: A.**
- **MI-Q25 ("every shipped set against the register").** A field of a non-silicon set counts as covered when one of
  these holds:
  - it has a 31.31 row;
  - it is 0 (not modelled, or the default form);
  - it is bit-equal to the legacy silicon value (inherited, and judged by that row).

  A frozen test asserts that every field the resolver knows is covered, for Ge, GaAs, InGaAs and 4H-SiC.
  AlGaAs's fields are GaAs's except those judged by MI-Q23.
  - **Ruling: A.** (B) would list the uncovered fields without asserting.
- **MI-Q26 (31.30's offset row).** Step 1's row "conduction_offset_eV" binds no shipped field: NiTCAD stores
  χ_ox = 1.05, which is silicon's χ less Gehring's offset.
  - (A) Rebind it as the χ_ox row. Its grade is L, because its input, silicon's χ 4.05, is L; Gehring's PA offset 3.0
    goes in the note.
  - Split the masses row into tunnel_mass_n and tunnel_mass_p, both tag C.
  - Add L rows for SiO₂'s κ 0.014 and ρc 1.628, whose inputs were not read.
  - (B) Leave the row unresolved (not_evaluated).
  - **Ruling: A.** 31.30 had not run.

Frozen expected outcomes (31.30, 31.31):
- No run error, fail or not_evaluated.
- Register rows:
  - pass only where the grade is PA and the value matches: SiO₂'s gap, and the radiative coefficients of InGaAs,
    4H-SiC and AlGaAs;
  - Ge's PA radiative row is a finding (a register mismatch, 6.4e-14 against 6.41e-14);
  - every other row is a finding (PC, L or C).
- Composition comparisons: findings (PC). Not_applicable where x lies beyond a law's printed range (N_c above 0.41).
- Coverage: no uncovered field.

**15.19 as built (2026-11-06).**

Files:
- `verify/reference/model_parameters.hpp`:
  - the 31.30 block revised by MI-Q26, now 7 rows;
  - 86 rows for 31.31;
  - the 7 AlGaAs `composition_laws`.
- `tvb31_common`:
  - `bound_field` resolves `radiative_cm3_s`;
  - `shipped_value` resolves the four constant sets, SiO₂ (`insulator_value`) and AlGaAs's B;
  - `non_silicon_sets()` and `field_names()` serve the coverage rule.
- `verify/benchmark/tvb31_materials.cpp`: `run_tvb31_30` and `run_tvb31_31`.
- `tests/verify/benchmark/tvb31_materials_test.cpp`.
- The registry holds 41 benchmarks.

Outcomes (Release, first run as frozen; under 0.1 s for the two):

| Entry | Pass | Finding | Not applicable | Not evaluated |
|---|---|---|---|---|
| 31.30 | 1 (the gap 9.0, PA) | 6 (ε_r, χ_ox, κ, ρc: L; both masses: C, pinned and matching) | 0 | 0 |
| 31.31 | 3 (B of InGaAs, 4H-SiC, AlGaAs: PA) | 152 | 1 (AlGaAs N_c at x = 0.45) | 0: every field is covered (MI-Q25) |

31.31's 152 findings break down as follows.

Printed values met: 21 PC rows.
- Ge: ε, β, and both μ bounds.
- GaAs: ε, E_g0, α, β, N_c, χ, both μ bounds, C_n (~1e-30) and B 7.2e-10.
- InGaAs: ε, N_c, N_v and both μ bounds.
- 4H-SiC: the hole μ bound 120 and κ 3.7.

Printed values missed: 26 rows, NiTCAD's value against the printed one.

| Material | Misses |
|---|---|
| Ge | E_g0 0.744 / 0.742; α 4.774e-4 / 4.8e-4 (NiTCAD's pair is Thurmond 1975's, not read); N_c 1.04e19 / 1.0e19; N_v 6.0e18 / 5.0e18; χ 4.13 / 4.0; T exponents −2.33 / −1.66 (electrons) and −2.23 / −2.33 (holes); Auger 2.8e-31 and 9.9e-32 / ~1e-30 (silicon's values carried over); P 0.012 / 0.013 eV; B 0.0104 / 0.011 eV; B_rad 6.4e-14 / 6.41e-14 (PA: a register mismatch in the coefficient the code cites Ioffe for) |
| GaAs | N_v 7.0e18 / 9.0e18; C_p 1e-31 / ~1e-30; Si 0.0058 / ~0.006 and C 0.0263 / ~0.02 eV (NiTCAD cites PRB 43, 14734, not read) |
| InGaAs | χ 4.55 / 4.5; Auger 1e-30 and 1e-31 / 7e-29 |
| 4H-SiC | α 3.3e-4 / 6.5e-4 and β 0 / 1300 K (NiTCAD's own comment calls its Varshni pair "tuned only to E_g(300 K)"; with β = 0 it gives E_g(300 K) = 3.201 eV against the printed 3.23); μ_n bound 950 / 900; τ_n0 5e-7 / ~1e-9 and τ_p0 5e-7 / ~6.0e-7 (the longest measured lifetimes); Auger 2.8e-31 and 9.9e-32 / 5.0e-31 and 2.0e-31 (silicon's values carried over) |

Other findings:
- **L, 36 rows:** no value read, findings by MI-Q5.
- **Composition laws, 69 points (PC):** 4 lie within the printed digits: ε, χ and μ_n at x = 0, and N_c at x = 0.15,
  where the two curves cross. The other 65 do not. At x = 0.45:

  | Quantity | NiTCAD | Ioffe's law |
  |---|---|---|
  | ε | 11.73 | 11.62 |
  | E_g(300 K) | 1.9836 | 1.9852 (1.5 meV lower at every x: the GaAs Varshni E_g(300 K) is 1.4225 against the printed 1.424) |
  | χ | 3.69 | 3.58 (the default share 0.68 against Ioffe's 0.88) |
  | N_v | 6.1e18 | 1.23e19 (NiTCAD's slope falls, Ioffe's rises) |
  | μ_n | 6025 | 125 |
  | μ_p | 332 | 83 |

  NiTCAD's own comment calls these scalings crude; their range, x ≤ 0.45, runs past Ioffe's N_c range x < 0.41 (the
  Γ–X crossover).

These are table findings (MI-Q10): they are reported, and none is fixed in TVR-3. The fix proposals go to the owner
together:
- Ge: the band set, Auger, and the dopant levels (two sources disagree: Ioffe against NiTCAD's uncited 12.0 / 10.4 meV).
- GaAs: N_v.
- InGaAs: χ and Auger.
- 4H-SiC: the Varshni pair, the lifetimes and Auger.
- AlGaAs: the composition laws.

Their effect on results is material-specific, and silicon is untouched.

Test tightening (after the first run, TVM2 allows it). A PC or L row is a finding whether it matches or not, so the
frozen counts could not see a changed value. The test now pins the split in the table above:
- 21 PC rows met and 25 missed;
- 1 PA mismatch;
- 36 L rows;
- 4 law points met.

L rows carry no value, so a change to an L-graded shipped value is invisible by design until a source is read.

Mutations (after a green baseline; filter `TVB-31.30*,TVB-31.31*` on `nitcad_verify_benchmark_test`). There were
11 one-literal changes to shipped values, each reverted by `git checkout`, and 10 were caught:
- M01: Ge ε 16.2 → 16.3 (PC met → missed);
- M02: GaAs N_v 7.0e18 → 9.0e18 (missed → met);
- M03: InGaAs B (PA pass → finding);
- M04: 4H-SiC κ;
- M05: the SiO₂ gap (PA);
- M06: the SiO₂ hole mass (tag C: the pin fails);
- M07: AlGaAs ε slope 2.6 → 2.84 (law points met);
- M09: Ge's hole N_ref (no row, ≠ silicon: coverage not_evaluated, MI-Q25);
- M10: AlGaAs B (PA).

**M08 survived:** AlGaAs N_c slope 0.5 → 0.6.
- No point changes between met and missed. The one met point, x = 0.15, stays within its ±1.6e16 tolerance.
- A law comparison is a finding (PC) either way, so it can only see a change that crosses the printed digits.
- This is recorded as the limit of a PC law check, not closed by a regression pin on NiTCAD's own numbers (that would
  be I0 evidence).

`git status src/physics include/NiTCAD/physics` is clean afterwards.

**15.20 Step 10: TVR-3 close-out (2026-11-06).**

**Debug cost and MI-Q9's split.** At step 9's end, the serial Debug verify (`ctest -R verify -LE full`, Debug,
serial) took 153.6 s. TVR-3's tests were about 91 s of that, against MI-Q9's allowance of 8 s and CV-Q19's 75 s
budget; 31.19 alone was 45 s. MI-Q9, ruled before code, applies: the dense grids move to `full`, and `fast` keeps
every switch point and one point a decade.
- `tvb31::per_decade(tier)` is 1 in fast and 64 in full and nightly.
- `tvb31::keep(k, n, tier)` keeps every 16th point and the last of a linear grid in fast. It is used in 31.05's
  explicit grid, 31.16's E_F, 31.17's η and depth sweeps, 31.22's d grid, 31.23's E, 31.24's x, 31.27's z, and 31.29's
  x grids.
- Four entries thin a set of non-switch points in fast:
  - 31.22's segment ends: the gap edges and their neighbours, not the interior ends 0.05–0.75. Its near-coincident
    pairs keep k = 3, 8 and 14, and its path z keeps 1e-10, 1 and 800.
  - 31.19's compensated offsets keep k = 1, 6 and 12.
  - 31.25 runs at 300 K and F₁ − F₂ = 0.1 and 1 eV.
- Every switch point, its `nextafter` neighbours, and every documented cut and range edge is added outside these grids,
  in every tier.
- Each fast point is therefore a point of the full grid, with the same comparison and outcome there.
- `nitcad verify` (default tier fast) runs the thinned grids. `--tier full` runs the frozen ones.

**MI-Q27 (where the frozen dense tests run).** Ruled A (owner standing instruction, "choose recommended",
2026-11-06).
- (A) As built:
  - The frozen tests of 31.06–31.31 call `Tier::full`. They carry the hidden Catch2 tag `[.tvb31_full]`, so the
    per-file ctests skip them.
  - The Release-only ctest `verify_benchmark.tvb31_full` runs them, with exactly their frozen counts. CI's release
    preset runs it in 11.7 s.
  - Both builds run `tvb31_fast_test.cpp`. It holds the frozen tests' per-comparison rules without their counts: no
    run error or fail, not_evaluated only by TB-Q4, notes only 31.25's diagnostics or a staleness note, and at least
    one comparison per entry but 31.25. A staleness note is allowed in fast, because a thinned grid may miss a
    finding's region; staleness is judged in full.
  - 31.01–31.05 take no tier, and their frozen tests run in both builds.
- (B) The `full` ctest label, which is local only (TVR-2's NITCAD_VERIFY_FULL_TIER). CI would never see the frozen
  TVB-31 outcomes.
- No criterion or expected outcome changed (TVM2): the full grids are the frozen ones.

Measured afterwards (serial Debug verify, the machine otherwise idle):
- 69.89 s overall, against CV-Q19's 75 s.
- TVR-3's share is 7.88 s, against MI-Q9's 8 s:

  | Tests | Time (s) |
  |---|---|
  | The seven L0 `model_*` tests | 2.08 |
  | `tvb31_band` | 1.36 |
  | `tvb31_fast` | 4.37 |
  | The seven per-file ctests, now empty | 0.07 |

**Suites:**
- Debug: `ctest -R verify -LE full`, 60/60 pass.
- Release, `ctest -j 6`, everything: 190 of 195 passed and none failed.
  - Claude Code stopped the run when the machine ran critically low on memory.
  - It was starting `verify_full.tvb04`, TVR-2's local-only full-tier ladder (about 20 GB). It runs only because
    this build tree sets NITCAD_VERIFY_FULL_TIER.
  - Not run: `verify_full.tvb04`, `app.run_controller`, `app_window.analysis_window`, `app_window.self_test` and
    `app_window_3d`. None touches TVR-3's code.
  - `verify_benchmark.tvb31_full`, `tvb31_fast`, and every other verify, physics, solve and deck test passed.

**TVB-31 outcomes (fast tier, `nitcad verify --id`), 31 entries:**
- pass, 9: 31.02, 31.03, 31.09, 31.13, 31.14, 31.16, 31.18, 31.23 and 31.26;
- finding, 19: 31.01, 31.04–31.08, 31.10–31.12, 31.17, 31.19–31.22, 31.24, 31.27, 31.28, 31.30 and 31.31;
- not evaluated, 2: 31.15 and 31.29 (TB-Q4's subnormal points, MI-Q19);
- unverified, 1: 31.25 (diagnostics only, MI-Q8).
- No entry fails.

**Provenance summary (R-3).** The grade of each entry is printed in the report's `reference` text; each register row's
is printed in its comparison's reason.
- Entries: PA 7 (31.02, 31.03, 31.18, 31.22–31.25); PB 9 (31.01, 31.06–31.09, 31.11, 31.20, 31.21, 31.27); PD 15 (the
  rest, which are derived here and checked against a second method or by the value rule on the register's numbers).
- Register (`model_parameters.hpp`), 193 rows: 29 PA, 41 PB, 62 PC, 55 L, 6 C; and 7 composition laws (PC).
- The L rows are the parameters with no source read. MI-Q6's primaries would raise them:
  - Caughey–Thomas and Baccarani–Ostoja (F3);
  - the lifetime set;
  - Dziewior–Schmid;
  - the Kane pair's origin (F5);
  - γ(T) (F4);
  - ε_r and χ;
  - the non-silicon saturation velocities and SiC's levels.

**Findings for the owner (MI-Q10: none fixed in TVR-3; each with a fix proposal in its section).** The rounding-floor
and cancellation findings, each a partial or value losing digits where its shipped form cancels:

| Finding | Where | Section |
|---|---|---|
| (vi) | n_i near 9 K | 15.13 |
| (vii) | the energy-balance μ₀ partials | 15.14 |
| (viii) | the SRH partials | 15.15 |
| (ix) | the trap-kinetics partials | 15.15 |
| (x) | the FD neutral η floor | 15.16 |
| (xi) | wkb_segment's short segments, and its midpoint truncation below D = 1e-10 | 15.17 |
| (xii) | the path's dG/dI_ik as z → 0 | 15.17 |
| (xiii) | dh/dx as x → −∞ | 15.18 |
| (xiv) | dN/dE where both levels lie well above E | 15.18 |

Their effects are absolute errors of order ε in quantities of order 1, or relative digits lost in partials near
degenerate inputs. No device result changes visibly.

Other findings:
- **Documented departures, the code against the published form (15.4, MI-Q4):**
  - F1: Slotboom's subtraction and cut.
  - F6: Masetti's T-scaled μ_max.
- **Attribution questions:** F2–F5 (15.4, 15.10).
- **Register mismatches:**
  - silicon: 15.10–15.18 (for example, the conductivity exponent 1.33 against the printed 1.3);
  - SiO₂ and the non-silicon sets: 15.19 (26 rows, and AlGaAs's laws).

**Next.** Unit 55's later phases (TVR-4/5) are not started. Every commit is local, on `verify/model-implementation`,
stacked on PR #78, which is open and not merged. Pushing, opening a PR or merging needs the owner's authorization.

#### 16. TVR-4 and TVR-5 (Units 55d, 55e) plan: assembled physics (proposed 2026-11-06; on `verify/device-theory`; rulings before code)

**Owner request (2026-11-06):** "plan for TVR-4 and TVR-5 (assembled physics)". TVR-3 checked every model *function*
(M class, 15.1). These two phases check the physics *as assembled and solved*: the device theories (A class), the
linear and transient responses, and the interface, generation and quantum terms. This section plans both phases. It
also reports, before any code, the places where §6's catalogue does not fit the code (16.2, 16.12). Each phase freezes
its criteria in a separate commit before the code that runs them (TVM2).

**Branch.** `verify/device-theory`, from `verify/model-implementation` at `dab9ec5` (PR #79, stacked on PR #78; neither
merged). TVR-4 needs TVR-2's first integral, KN reference and study driver, and TVR-3's re-coded models. It is rebased
on `main` once #78 and #79 merge. TVR-5 gets its own branch, stacked on this one, when TVR-4 closes (16.14). Local
only.

##### TVR-4 (Unit 55d): device theory, small signal, transient, thermal

**16.1 Scope.** In:
- (a) The benchmarks of §9's TVR-4 row (TVB-10, 13, 15–21, 27, 28), plus TVB-05 and 06, which CV-Q3 B deferred from
  TVR-2.
- (b) The audit's TVR-4 items (§10.7): R-4's reverse-bias dQ/dV, R-5 (the diode family with consistent assumptions),
  R-9 (the Shockley admittance with both carriers and C_T, eq 4.23) and R-10 (Pao–Sah and the charge sheet).
- (c) The new references, base-only (TVD-19), each with its L0 tests and its contract (§3.3).
- (d) The two changes to existing code the catalogue needs (16.2 (a), (b)), each only if ruled, each with its own
  tests, each bit-identical for every existing device and run.
- (e) The TVB-32 consistency suite on every new benchmark device (§6: "every P1/P2 benchmark device").
- (f) Every benchmark registered in TVR-1's registry with its tier, and its frozen-outcome test.

Out:
- TVR-5's benchmarks (16.12).
- The existing tests that TVB-16, 17, 18, 20, 21, 27 and 28 overlap (`bias:165`, `:284`, `small_signal:856`, `:994`,
  `mosfet:222`, the electrothermal U4/DC7 gates). They stay as they are (TVD-14). "Re-homing" (§6) means a new,
  independent benchmark beside the original; the original's gate is kept.
- Any change to NiTCAD's physics, tables, tests, gates or tolerances. A miss is a finding, and any fix is proposed
  separately (as Unit 56 did).

**16.2 Found before code (the catalogue against the code).**
- **(a) Insulator-only devices are refused** (`device.cpp:75`, "device needs at least one semiconductor region"), and
  the reference node assumes a semiconductor node (`device.cpp:472`). TVB-05 and 06 cannot be built without a
  device change. CV-Q3 deferred them here, so they must now be built (DT-Q2 A) or recorded as unbuildable (DT-Q2 B).
- **(b) The transient API cannot start from an excess-carrier state.** `transient.hpp:83`: "`initial` is a starting
  guess for the steady state at t = 0". TVD-15 said the initial state "exists"; it does not. Nor does prescribed
  generation (§7). So TVB-19, dielectric relaxation and the Haynes–Shockley pulse have no way to inject carriers.
  Options in DT-Q3.
- **(c) TVB-17's "exact dQ/dV at reverse bias" is exact only for a model NiTCAD does not solve.** TVR-2's first
  integral is exact at equilibrium, where the quasi-Fermi levels are flat. Under bias the first integral holds for the
  *quasi-equilibrium* (QE) model: electrons referred to φ_n of the n side, holes to φ_p of the p side, both flat
  across the junction. Drift-diffusion departs from QE by the quasi-Fermi drop the current causes. That departure is
  small where the carriers' charge matters and large only where their density is negligible. TVB-17 is therefore E
  against the QE model, with a predicate on the solved state (DT-Q4).
- **(d) The diode laws need consistent assumptions (I-2, R-5).** Shockley's law is the whole current only without
  depletion-region recombination, and its long-base form needs recombination in the neutral regions. So TVB-16 has
  two configurations:
  - (i) SRH, Auger and radiative off: only the short-base law applies (the contacts are the only sinks). W(V) comes
    from the QE first integral (not W(0)), and the base length is measured from the depletion edge.
  - (ii) SRH on: the long- and short-base laws with L = √(D τ), plus the Sah–Noyce–Shockley (SNS) integral of U(ψ)
    over the *reference's* QE potential (not NiTCAD's field: C-2's lesson).
  - The reverse-leakage exponent (I-3) is judged in (ii), against the SNS integral, not a fixed band.
  - Mobility is constant (uniform doping per side, field mobility off) and its value and the lifetimes are given in
    the spec (R-7), so D = V_T μ and τ enter the reference without NiTCAD's tables. NiTCAD's tables must agree with
    the spec values; a mismatch is an L1 finding (§2.2), not hidden by sharing.
- **(e) TVB-15's injection.** A minority-diffusion region with SRH needs carriers injected. With no prescribed
  generation (b), a junction would inject them, and the boundary condition would then be the law of the junction (A,
  not E). The TVQ3 hook already supplies a *static*, state-independent generation G(x). With a localized G in a long
  uniform p-type bar, the linear minority equation D Δn'' − Δn/τ = −G has a closed-form solution (Green's function,
  cosh/sinh), with τ the *exact* linearization of SRH at the equilibrium state (it includes the n₁, p₁ terms),
  computed in the reference from the spec. NiTCAD solves full DD, not the linear model, and the two differ by
  O(Δn/N_A) (the ambipolar field and the SRH nonlinearity). Hence two comparisons:
  - E, with the predicate max Δn/N_A ≤ 1e-8 measured on the solved state: the model error is then about 1e-8, two
    decades below the 1e-4 tolerance and below the finest discretization error, so it cannot pollute p̂.
  - A, at Δn/N_A from 1e-3 to 1e-1, with ε_approx derived in step 1 and frozen before it runs.
- **(f) TVB-13 is exact in a different way than the catalogue says.** In a uniform resistor with ohmic contacts,
  n = N everywhere, so J = q N μ(E) E forces E to be uniform. Canali drift is then a D-class identity, as TVB-11 is.
  With energy balance, T_n = T_L at the contacts, so there is a boundary layer of length λ ≈ v τ_w (an upper bound on
  λ derived in step 1, PD). In the interior the uniform state solves the discrete equations exactly, so there is no
  discretization error to converge, only the boundary layer's tail: T_n is judged as D, not E, on the nodes whose
  distance d from both contacts gives e^{−d/λ} ≤ 1e-13 (a predicate with its measured value), and no order is
  computed.
- **(g) TVB-10's frequencies and its HF model.** HF C-V needs f above the minority carriers' generation response and
  below their redistribution across the inversion layer. The HF reference is therefore not "minority density frozen"
  but the redistribution limit: the minority quasi-Fermi level spatially uniform, its total number conserved, the
  majority carriers at equilibrium (to be confirmed against the source in step 1). NiTCAD's `frozen_minority_
  capacitance` is not used (C-4: I0). The plateau is measured and must be flat to the tolerance across a decade
  (else NOT APPLICABLE). Deep depletion needs a transient ramp faster than generation. Its reference is the KN charge
  with zero minority carriers, exact for that state, under the predicate minority charge ≤ 1e-3 |Q_s| on the solved
  state.
- **(h) TVB-21 needs no separate formula.** S = ln10 V_T (1 + C_dep/C_ox) is the depletion approximation *of* the
  Pao–Sah current in weak inversion. The A reference is S from the Pao–Sah reference itself (TVB-20) at the same
  V_G. The closed formula is reported beside it as a diagnostic.
- **(i) TVB-20's assumptions.** Pao–Sah requires constant mobility (field mobility and Lombardi off), uniform
  substrate, Boltzmann, no BGN, a lumped gate, and the gradual channel. Gradual channel is checked on the solved
  state: max |∂E_x/∂x| / |∂E_y/∂y| in the channel ≤ a frozen bound, with source and drain far enough apart.
- **(j) TVB-05's spherical capacitor cannot be meshed.** `Geometry3D` is a convex polyhedron with planar interfaces
  (`generate.hpp:263`). A spherical shell has a curved, interior boundary, so it cannot be built. The coaxial
  capacitor can be built in 2D (per unit length): the 2D generator takes region polygons that tile the domain without
  nesting, so the annulus is tiled by sectors of one region, with the inner and outer circles faceted. The spherical
  case is recorded as unbuildable (unverified, not failed).
- **(k) TVR-0's C-9 for thermal.** The existing electrothermal gates take expected values partly from NiTCAD's
  physics functions. The re-homed TVB-27 and 28 take κ(T) = κ₃₀₀ (T/300)^−a (`thermal.hpp:4`'s form), S and Π from
  the benchmark spec, re-coded (R-7).
- **(l) TVB-18 needs recombination, which brings a term eq 4.23 omits.** The long-base admittance needs SRH on, and
  SRH also recombines in the depletion region. Its admittance share is bounded in ε_approx (step 1), and the forward
  biases judged are those where the reference's SNS current is ≤ 1e-2 of the diffusion current (a predicate).
- **(m) TVB-19's linear models include diffusion.** The excess majority charge spreads while it relaxes, but its
  integral decays exactly as exp(−t/τ_d) in the linear model (ohmic ends far from the pulse, their flux bounded in
  the predicate). The decay rate is judged on that integral. The Haynes–Shockley pulse is E only under the same
  low-injection predicate as TVB-15 (Δp/N_D ≤ 1e-8).
- **(n) Pao–Sah's inner integral.** The behaviour of n/E at the bulk end of the inner ψ integral (where E → 0) is
  derived in step 1 before the quadrature is chosen. The quadrature does not hide a divergence.
- **(o) Reach of the two code changes.** (a) and (b) may reach beyond `device.cpp` and the hook: assembly with no
  carrier rows, the equilibrium and bias solves, and the transient passing t to the residual (the hook's profile adds
  its own breakpoints to the step-landing set). Each step measures the reach first and stops for the owner if the
  change touches more than the device check, the reference-node choice, the hook and the transient's breakpoint set.
- **(p) Error interactions (TVD-17, R-12).** Every E and A benchmark runs the 2×2 (h, tol) factorial at its finest
  level and reports ε_int, and the transients add Δt (2×2×2). An ε_int above 10% of U blocks an E pass (§4.1).

**16.3 New references (`verify/reference`, base only; L0 tests in each one's commit).**

| Ref | Contents | Method; second method | Used by |
|---|---|---|---|
| capacitor | C = 2πεL/ln(b/a); series stack C⁻¹ = Σ t_i/ε_i, D continuous | Closed form; the faceted polygon's own capacitance by Richardson extrapolation (the O(h²) faceting term is measured, not assumed) | TVB-05, 06 |
| minority_diffusion | D n'' − n/τ = −G on [0, L] with ohmic ends (n = n₀); G a Gaussian or box | Green's function, integrated by Gauss–Kronrod in dd; Chebyshev collocation | TVB-15 |
| qe_junction | TVR-2's first integral under bias in the QE model: ψ(x), the edges, Q(V) and C = dQ/dV | tanh–sinh and Brent (as TVB-07); C by complex step on Q(V); Chebyshev | TVB-16, 17, 23 (TVR-5) |
| diode | Shockley J₀ (long and short base, the base from the QE edge); the SNS integral of U(ψ) over the QE profile | Closed form and Gauss–Kronrod; the SNS integral also by tanh–sinh | TVB-16 |
| admittance | Shockley eq 4.23: Y = Σ_{n,p} G₀ √(1 + iωτ) (long base) or the coth form (short), plus iω C_T | Closed form in complex dd; the linearized minority ODE solved directly (Chebyshev, complex) | TVB-18 |
| mos_hf | KN charge in the HF redistribution limit (16.2 (g)), and with no minority carriers (deep depletion) | Closed form and Brent; HF C by complex step | TVB-10 |
| pao_sah | I_D = q μ (W/L) ∫_{V_S}^{V_D} ∫ n(ψ, V)/E(ψ, V) dψ dV, the inner integral from the bulk to ψ_s(V_G, V), ψ_s by Brent on TVB-09's KN; Brews' charge sheet; S | Nested adaptive Gauss–Kronrod; tanh–sinh at the endpoint whose behaviour step 1 derives (16.2 (n)) | TVB-20, 21 |
| linear_transient | Dielectric relaxation: the excess charge's integral as exp(−t/τ_d), τ_d = ε/(q μ N); Haynes–Shockley: a drifting, spreading, decaying Gaussian convolved with the source's time profile | Closed form; the convolution by Gauss–Kronrod in dd; the PDE by method of lines with Dormand–Prince (TC6) | TVB-19 |
| thermal_slab | Kirchhoff transform for κ = κ₀ (T/T₀)^−α; Joule parabola; Seebeck and Peltier integrals | Closed form; the ODE by Dormand–Prince | TVB-27, 28 |

The homogeneous drift and energy-balance references (TVB-13) are TVR-3's re-coded `model_mobility` and `model_thermal`,
reused unchanged.

**16.4 The benchmarks.**

| ID | Construction | Reference; class | Judged | Families (fast / full) |
|---|---|---|---|---|
| TVB-05 | Coaxial capacitor, 2D generated annulus (sectors), electrodes on both circles (DT-Q2 A) | capacitor; E | C per unit length, p̂ (faceting included in p_th 2) | — / 4 levels |
| TVB-06 | SiO₂ / high-k stack between electrodes, 1D and 2D tensor, interfaces on dual faces | capacitor; D | C and the field ratio ε₁/ε₂ | 1D N = 32 / + 2D |
| TVB-10 | MOS-C 1D, SRH on: HF C-V by small signal at f in the measured plateau; deep depletion by a transient V_G ramp | mos_hf; A (HF), E (deep depletion, with the predicate) | C_HF(V_G), C_min; Q_s in deep depletion | — / 1D ladder |
| TVB-13 | Uniform resistor, field mobility on (Canali); then energy balance on | model_mobility, model_thermal; D (drift, and T_n on the judged interior, 16.2 (f)) | v = J/(qN) against v(V/L); T_n against the homogeneous T_n(E) | 1D (fast: drift) / + energy balance |
| TVB-15 | Long uniform p-bar, ohmic ends, SRH only, the TVQ3 hook's localized G | minority_diffusion; E (Δn/N_A ≤ 1e-8) and A (Δn/N_A 1e-3 to 1e-1) | relative Δn in L∞; p̂ | — / 1D ladder |
| TVB-16 | Abrupt 1D diode, constant μ; (i) no G−R, short base; (ii) SRH on, long and short base | diode; A | J(V) by log error; ideality; the reverse exponent (ii) | — / 1D ladder per V |
| TVB-17 | Abrupt 1D diode at reverse bias, no G−R | qe_junction; E (QE, predicate DT-Q4); the ε/W form's difference reported | C_j(V), p̂ | — / 1D ladder |
| TVB-18 | Long and short diode, forward bias, low injection | admittance; A | Y(ω) magnitude and phase over 1e-3 ≤ ωτ ≤ 1e3 | — / 1D |
| TVB-19 | Uniform n-bar: (i) a majority-charge pulse (dielectric relaxation); (ii) a minority pulse in a uniform field (Haynes–Shockley); injection per DT-Q3 | linear_transient; E | decay rate; pulse centre, width and area; time-step p̂ (BDF2: 2; BE: 1) and the 2×2×2 factorial (h, tol, Δt) | — / 1D |
| TVB-20 | Long-channel n-MOSFET, 2D tensor, constant μ, lumped gate | pao_sah; A (gradual channel); Brews as second A | I_D(V_G, V_D), V_T, g_m | — / 2D ladder |
| TVB-21 | TVB-20's device in weak inversion, low V_D | pao_sah; A | S | — / with TVB-20 |
| TVB-27 | 1D slab with κ(T) from the spec, isothermal ends; Joule-heated resistor | thermal_slab; E | T(x) in L∞, p̂ | — / 1D ladder |
| TVB-28 | Seebeck and Peltier: a 1D bar with ∇T; the energy identity | thermal_slab; E / D | S, Π; the energy residual | — / 1D |

Every reference input comes from the benchmark spec, never from NiTCAD's tables (§2.2, R-7). Every comparison carries
its applicability predicates (TVD-6), and a false one is NOT APPLICABLE with the predicate and its measured value
named (TVD-16).

**16.5 Sources (step 1, read before criteria freeze; grades per §3.3).**

| Source | Status | Needed for |
|---|---|---|
| Shockley, BSTJ 28 (1949), eqs 4.23–4.24 | Read in TVR-0 (OCR); PA | TVB-16, 18 |
| Sah, Noyce and Shockley, Proc. IRE 45 (1957) | To read | TVB-16 (SNS) |
| Pao and Sah, Solid-State Electronics 9 (1966) | To read (paywalled per §9); PB via textbook reproductions otherwise | TVB-20 |
| Brews, Solid-State Electronics 21 (1978) | To read; as above | TVB-20, 21 |
| Haynes and Shockley, the drift-mobility experiment (Phys. Rev., 1951) | To locate and read | TVB-19 |
| The HF MOS C-V derivation (Nicollian and Brews, *MOS Physics and Technology*, ch. 4, or an open reproduction) | To locate | TVB-10 |
| Kingston and Neustadter (TVR-2) | As graded in TVR-2 | TVB-10, 20 |
| Wachutka, IEEE Trans. CAD 9 (1990) (thermodynamic heat generation) | To read | TVB-28 |
| Kirchhoff transform (Carslaw and Jaeger, or any heat-conduction text) | Derived here (PD) with a second method | TVB-27 |
| Canali et al. (1975) | Read in TVR-3 (15.10) | TVB-13 |

Where only PB is reachable, the grade is printed in every report line. A PC-only source makes the comparison
"provenance insufficient" and it cannot pass (TVD-18).

**16.6 Order of work (one commit each; L0 tests before the benchmark that uses them; mutation runs per commit as in
TVR-3).**
1. Reading log and the derivations (PD, each with its second method): the bound on λ (TVB-13); the HF
   redistribution model (TVB-10); the behaviour of Pao–Sah's inner integrand at the bulk end; ε_approx for TVB-10
   HF, 15, 16, 18, 20 and 21. Recorded as their own subsection; then the A-class criteria freeze (DT-Q8).
2. The insulator-only device change (DT-Q2 A) with its tests; capacitor reference; TVB-06, then TVB-05.
3. The hook's time profile (DT-Q3 A) with its tests.
4. minority_diffusion; TVB-15.
5. qe_junction; TVB-17.
6. diode; TVB-16.
7. admittance; TVB-18.
8. mos_hf; TVB-10.
9. pao_sah; TVB-20, TVB-21.
10. linear_transient; TVB-19.
11. TVB-13.
12. thermal_slab; TVB-27, TVB-28.
13. TVB-32 consistency (C1–C12) on every new device; findings to the owner.
14. Close-out: timing, suites, provenance summary, ARCHITECTURE 10 catalogue links (documentation only).

**16.7 Tiers and time.** Serial Debug `nitcad verify` is at 69.89 s of CV-Q19's 75 s (15.20). TVR-4 adds to `fast` only
its two D-class cases (TVB-06 1D, TVB-13 drift), measured, with ≤ 2 s total. Everything else is `full`. The frozen full
outcomes run in a Release-only ctest `verify_benchmark.tvr4_full` (as `tvb31_full`, so CI sees them) if it takes ≤ 120 s
on CI. Ladders beyond that budget go under `NITCAD_VERIFY_FULL_TIER` (local, as TVB-04's).

**16.8 Criteria for the D and E comparisons (frozen with these rulings, 2026-11-06, before code; TVM2).** Notation as
14.7: ε = 2⁻⁵², N the nodes along the path, p̂ by §4.2 on r = 2 ladders, the CV-Q8 iterative rule, TB-Q4's
reference-uncertainty rule.

| ID | D criterion | E (a) order | E (b) accuracy on the finest level | Ladder |
|---|---|---|---|---|
| TVB-06 | C and the field ratio, relative ≤ 64 ε N² | — | — | 1D N = 32, 128; 2D 32² |
| TVB-05 | — | C: p̂ ∈ [1.8, 2.3] | C relative ≤ 1e-3 | 4 levels, 2D generated |
| TVB-13 drift | v relative ≤ 64 ε N² at V/L from 1e2 to 1e5 V/cm | — | — | 1D N = 64 |
| TVB-13 T_n | T_n − T_L relative ≤ 64 ε N² on the judged nodes (16.2 (f)) | — | — | 1D N = 128, 512 |
| TVB-15 linear | — | relative Δn in L∞: p̂ ∈ [1.8, 2.3] | ≤ 1e-4 | 1D N = 64…512 |
| TVB-17 | — | C_j: p̂ ∈ [1.8, 2.3] | relative ≤ 1e-3 | 1D N = 64…512, at each V |
| TVB-19 | — | in Δt at fixed h: p̂ ∈ [1.8, 2.3] (BDF2), [0.8, 1.3] (BE); in h at fixed Δt: [1.8, 2.3] | decay rate relative ≤ 1e-3; pulse centre ≤ 1e-3 of the drift length, width and area relative ≤ 1e-3 | Δt and h ladders, 4 levels each; the 2×2×2 factorial at the finest |
| TVB-27 | — | T in L∞: p̂ ∈ [1.8, 2.3] | relative ≤ 1e-4 | 1D N = 64…512 |
| TVB-28 | energy residual within C5's frozen bound (11.6), with C5's accepted Newton-term limitation printed | S, Π: p̂ ∈ [1.8, 2.3] | relative ≤ 1e-4 | 1D N = 64…512 |
| TVB-10 deep depletion | — | Q_s: p̂ ∈ [1.8, 2.3] | relative ≤ 1e-3 | 1D N = 64…512 |

A-class criteria (TVB-10 HF, 15 DD-vs-linear, 16, 18, 20, 21) are §4.4's A row, |q_sim − q_ref| ≤ U_num + U_ref +
ε_approx, with each ε_approx frozen in step 1 from its derivation and source, before the code that runs it (DT-Q8).

**16.9 Owner questions (recommendations in the last column).**

| # | Question | Recommendation |
|---|---|---|
| DT-Q1 | Scope as 16.1, including TVB-05/06 and R-4, R-5, R-9, R-10; the overlapping existing tests kept as they are? | Yes |
| DT-Q2 | Insulator-only devices: (A) a narrow change to `device.cpp:75` and `:472`: allow a device with no semiconductor region if every connected component of its mesh has an electrode, with the reference node an electrode node; its own tests; bit-identical for every existing device. (B) Record TVB-05/06 as unbuildable (unverified). | A. CV-Q3 deferred them here on the understanding that they would be built; B would leave the Poisson operator's only exact closed-form capacitances unchecked |
| DT-Q3 | Excess-carrier injection for TVB-19: (A) give the TVQ3 hook an optional time profile w(t), piecewise constant on the transient's waveform breakpoints (so steps land on its jumps); null by default, verification-only, not in the deck, bit-identical when absent. The pulse is then a short generation burst, which is what the experiment does physically, and the reference convolves with w(t). (B) A transient option to take `initial` as the t = 0 state instead of a guess: an API change to a user-facing function. (C) Defer TVB-19 | A. It extends a hook that is already verification-only; B changes user-facing semantics |
| DT-Q4 | TVB-17: E against the QE model, applicable only where the solved state's quasi-Fermi drop across the region in which the carrier's charge exceeds 1e-6 of the local doping is ≤ 1e-3 V_T (else NOT APPLICABLE with the measured drop)? | Yes |
| DT-Q5 | TVB-16 in two configurations (16.2 (d)), with constant μ from the spec; the 1.321e-2 short-base value of `bias:144` cross-checked, not gated? | Yes |
| DT-Q6 | TVB-21: S from the Pao–Sah reference (A), the closed formula a diagnostic? | Yes |
| DT-Q7 | TVB-05: coaxial only; the spherical capacitor recorded as unbuildable (no curved 3D geometry)? | Yes |
| DT-Q8 | Criteria: freeze 16.8 (D and E) with these rulings; freeze each A-class ε_approx in step 1, after its source is read and its derivation written, before the code that runs it? | Yes |
| DT-Q9 | Tiers as 16.7 (two D cases in `fast`, ≤ 2 s; `verify_benchmark.tvr4_full` Release-only if ≤ 120 s on CI; heavier ladders local)? | Yes |
| DT-Q10 | A failed predicate is NOT APPLICABLE; a benchmark whose every comparison is NOT APPLICABLE is "unverified" (TVD-16), never adjusted to apply? | Yes |
| DT-Q11 | Findings: none fixed in TVR-4; each gets a fix proposal and goes to the owner (MI-Q10's rule)? | Yes |
| DT-Q12 | Unit numbering and branch: Unit 55d on `verify/device-theory`; TVR-5 as Unit 55e on its own branch after TVR-4 closes? | Yes |

##### TVR-5 (Unit 55e): interfaces, generation, tunnelling, quantum

**16.10 Scope.** In: §9's TVR-5 row: TVB-08, 22–26, 29, 30, 33; the C10 limits left for TVR-5 (MLDA → classical far
from the interface, Lombardi → bulk as E⊥ → 0; §11.6); the new numerics TVD-22 placed here (the Sturm-sequence
eigensolver); TVB-32 on every new device. Out: as 16.1.

**16.11 The benchmarks (each re-planned in TVR-5's own proposal, which freezes its criteria).**

| ID | Construction | Reference; class |
|---|---|---|
| TVB-08 | Equilibrium abrupt heterojunction with an ε step, both interface modes (continuous DD and thermionic emission; their equilibria must agree) | First integrals on each side matched by D continuity (Anderson rule); E. `heterojunction:147` kept |
| TVB-22 | MOS-C with one trap level: G_p/ω and C_it(ω) | Single-level equivalent circuit, derived here (PD) from the small-signal SRH trap model, since eq 17 of Nicollian–Goetzberger was illegible (I-10); E for the model |
| TVB-23 | Abrupt 1D diode, impact ionization on | The ionization integral of TVR-3's re-coded vOdM over the *reference* field (qe_junction): E for M(V) below breakdown, under a low-current predicate (the carriers' space charge negligible, measured); A for NiTCAD's V_BR on Unit 19's traced branch against the reference's integral = 1 |
| TVB-24 | p-i-n at reverse bias (field near uniform, predicate measured); local Kane, then nonlocal paths | ∫ G(E(x)) dx on the reference field; the path integral by quadrature on the reference potential; E for the model as assembled |
| TVB-25 | MOS-C with gate tunnelling (Unit 51) | (1) E: Tsu–Esaki with WKB on the reference's own band diagram (KN): the assembled coupling; (2) H: WKB vs exact transmission (Airy functions, Gundlach) by oxide-field regime, recorded only (TVQ5); (3) A: Fowler–Nordheim and Simmons in their regimes |
| TVB-26 | Schottky diode (Unit 48) | Thermionic emission and the Crowell–Sze thermionic-diffusion combination; A. Image-force lowering NOT APPLICABLE (16.12 (a)) |
| TVB-29 | Cylindrical pn junction: 2D generated, a faceted circular junction | Radial Poisson–Boltzmann by Chebyshev collocation; E. Spherical unbuildable (16.12 (b)) |
| TVB-30 | MLDA vs 1D Schrödinger–Poisson, infinite and finite barrier | In-reference Sturm bisection + inverse iteration in dd (TVQ4 ruling), Eigen secondary; H, recorded only |
| TVB-33 | Aste–Vahldieck ballistic diode | P; see 16.12 (c) |

**16.12 Found before code (TVR-5).**
- **(a) Image-force lowering and Schottky tunnelling are not modelled** (`contact.hpp:22`). TVB-26's barrier-lowering
  comparison is NOT APPLICABLE with that predicate named. Unit 48's five accepted limitations are printed beside the
  outcome.
- **(b) No curved 3D geometry** (16.2 (j)). TVB-29's spherical junction is unbuildable; the cylindrical junction is
  built in 2D, its faceting measured as in TVB-05.
- **(c) TVB-33 compares different models (I-8).** NiTCAD's τ_w is constant (`energy_balance.hpp:26`); the paper's
  depends on T (its eqs 15–16), and its β is not stated. A P-class comparison would judge the model mismatch, not the
  code. The comparison is NOT APPLICABLE ("the model matches the publication's": false), with NiTCAD's peak velocity
  printed as a diagnostic. It becomes applicable only if a T-dependent τ_w is added, which is not in scope.
- **(d) New numerics.** Airy functions Ai, Bi and their derivatives in dd, scaled for large arguments (Bi grows as
  exp(2ζ^{3/2}/3)), with L0 against DLMF §9's tables and identities (the Wronskian); the exact transmission is
  cross-checked by integrating the Schrödinger equation through the barrier with Dormand–Prince (I3, two methods).
  The Schrödinger–Poisson reference's self-consistency loop must give its own U_ref (a residual bound, not an
  iteration count).
- **(e) C-2's circularity is avoided by construction.** TVB-23 and 24 integrate on the *reference* field, never on
  NiTCAD's. Findings (xi) and (xii) of TVR-3 are carried into TVB-24's report.
- **(f) Size.** TVR-5 is XL (§9). Split it into two units, each with its own proposal, rulings and criteria freeze:
  TVR-5a (TVB-08, 22, 23, 24, 26, 29; classical interfaces and generation) and TVR-5b (TVB-25, 30, 33 and the C10 MLDA
  limit; tunnelling and quantum). The Lombardi limit goes with TVR-5a's MOS devices.
- **(g) Sources, all to read in each sub-unit's step 1:** Anderson (Solid-State Electronics 5, 1962); Nicollian and
  Brews (book) for the trap circuit; Sze and Gibbons (1966) and Unit 19's records for breakdown; Kane's
  band-to-band rate (attribution still open: F5 of 15.4) and Unit 51's records; Gundlach (Solid-State Electronics 9, 1966); Simmons (J. Appl. Phys. 34, 1963); Crowell and Sze
  (Solid-State Electronics 9, 1966); Paasch and Übensee (MLDA); Stern's self-consistent inversion-layer papers as the
  SP cross-check. Grades as §3.3; a PC-only reference cannot pass.

**16.13 Owner questions for TVR-5 (the plan only; each sub-unit's proposal asks its detailed questions).**

| # | Question | Recommendation |
|---|---|---|
| AP-Q1 | Split TVR-5 into 5a (classical interfaces and generation) and 5b (tunnelling and quantum), as 16.12 (f)? | Yes |
| AP-Q2 | TVB-26's image-force and TVB-29's spherical case recorded as NOT APPLICABLE / unbuildable rather than adding the models or geometry? | Yes. Adding them is new physics and meshing, outside a verification unit |
| AP-Q3 | TVB-33 NOT APPLICABLE with the mismatch named, NiTCAD's peak printed as a diagnostic (16.12 (c))? | Yes |
| AP-Q4 | TVB-25 in three layers: E for the assembled WKB coupling, H (record only) for WKB vs Airy, A for FN and Simmons? | Yes |
| AP-Q5 | TVB-23/24 integrate on the reference field only (16.12 (e))? | Yes |
| AP-Q6 | Criteria for TVR-5 frozen in each sub-unit's proposal, after its step 1 reading, before its code? | Yes |

**16.14 Sequence.** TVR-4 steps 1–14 (16.6), close-out, the owner's review; then TVR-5a's proposal (rulings, criteria
freeze, code), then TVR-5b's. Nothing is pushed without the owner's authorization.

**16.15 Review of the plan (owner: "choose recommended and review the plan", 2026-11-06).** The first draft was
checked against the code and against §4's rules before the rulings. The review changed:
- TVB-16 (i): with no G−R there is no long base. Only the short-base law applies; the long base and the reverse
  exponent moved to (ii).
- TVB-15: "E for the linear model" cannot be judged on a full-DD solve unless the model error is far below the
  discretization error. E now has the predicate Δn/N_A ≤ 1e-8, and the A comparison runs at higher injection. τ is
  SRH's exact linearization, not τ_n.
- TVB-13: the interior energy-balance state is exact in the discrete equations, so an order study measured only the
  boundary layer's tail. It is D on the judged nodes, with a tail predicate.
- TVB-10: "minority density frozen" is not the HF limit of a 1D MOS-C small-signal solve. The reference is the
  redistribution limit (uniform minority quasi-Fermi level, number conserved), to be confirmed in step 1.
- TVB-18: SRH in the depletion region adds admittance that eq 4.23 omits. It is now bounded by a predicate and in
  ε_approx.
- TVB-19: the dielectric-relaxation decay is judged on the excess charge's integral (diffusion conserves it). The
  pulse needs TVB-15's low-injection predicate.
- TVB-20: the first draft put the inner integral's singularity at the surface. Where n/E is singular (at the bulk
  end, where E → 0) is derived in step 1 before the quadrature is chosen.
- TVB-28: the energy residual uses C5's frozen bound, not a bare 64ε, since C5's Newton-term limitation applies.
- TVB-23: E for M(V) needs a low-current predicate. V_BR is NiTCAD's traced value against the reference integral.
- Added: the 2×2 (h, tol) factorial on every E/A benchmark (TVD-17, R-12); the reach of the two code changes measured
  before they are made, with a stop for the owner (16.2 (o)); registry entries (16.1 (f)).
- Removed: a Kane date the register does not support (F5 is open).

Not changed, and checked:
- `device.cpp:75` and `transient.hpp:83` read as quoted.
- `Geometry3D` is convex with planar interfaces.
- `contact.hpp:22` states that there is no image-force lowering or tunnelling.
- The hook already has separate electron and hole generation (`models.hpp:23`), so a majority-only pulse needs no
  new source kind.
- The fast-tier budget: 69.89 s + ≤ 2 s < 75 s.

Open risks, for the owner:
- (1) DT-Q2 A may reach the solves (16.2 (o)).
- (2) Pao–Sah, Brews and Sah–Noyce–Shockley may be reachable only at PB.
- (3) TVB-10 HF's redistribution model is this plan's derivation until the source is read. If the source differs,
  the source wins and the change is recorded.

**16.16 Rulings (owner, "choose recommended", 2026-11-06).** DT-Q1 through DT-Q12 and AP-Q1 through AP-Q6 are
approved as recommended, on the plan as reviewed in 16.15:
- TVB-05 (coaxial only) and TVB-06 are built in TVR-4 with the narrow insulator-only device change (DT-Q2 A).
- TVB-19's injection is a time profile on the TVQ3 hook (DT-Q3 A).
- TVB-17 is E against the QE model under DT-Q4's predicate.
- 16.8's D and E criteria are frozen now. Each A-class ε_approx is frozen in step 1, before the code that runs it.
- TVR-5 is split into 5a and 5b. Each gets its own proposal and criteria freeze.
- No code was written with this plan. Step 1 (reading and derivations) is next, on the owner's "go".

**16.17 Step 1: reading log (owner "go", 2026-11-06).** Open copies were fetched with curl and rendered page by
page (WinRT PDF render). Nothing was taken from memory without a source line below.

| Source | How read | Grade | What it settles |
|---|---|---|---|
| Haynes and Shockley, Phys. Rev. 81, 835 (1951) | Scan pp. 835–838 (ibm1401.computerhistory.org) | PA (model statement); formula PD | §III: with the emitter current small enough that conductivity modulation is neglected, the pulse fronts have "the customary integral of the error function" shape. Diffusion spreads them equally ahead and behind. Footnote 9: recombination multiplies by exp(−t/τ). §VI: N = N₀ exp(−t/τ). This is the linear drift–diffusion–decay equation under low injection. Its Gaussian Green's function is derived here (PD), and conductivity modulation (Herring, BSTJ 28, 401) is the departure that the low-injection predicate excludes. |
| Pao and Sah, Solid-State Electronics 9, 927 (1966) | Primary paywalled. (1) H.-S. P. Wong, Stanford EE 316 slides 3-8 to 3-11 (hosted at dunham.ece.uw.edu); (2) Galup-Montoro, Schneider and Cunha, NSTI-Nanotech 2007, 474, eqs 1–2 | PB | Both give I_D = (μW/L)∫_{V_S}^{V_D}(−Q_I(V))dV. (1) gives Q_I = −q∫_{ψ_B}^{ψ_s}(n_i²/N_a)e^{(ψ−V)/V_T}/E(ψ,V)dψ, **with the inner lower limit ψ_B (the bulk Fermi potential), not 0**, and E² = (2kTN_a/ε)[(e^{−ψ/V_T} + ψ/V_T − 1) + (n_i²/N_a²)(e^{−V/V_T}(e^{ψ/V_T} − 1) − ψ/V_T)]. (2) uses the same double integral. Slide 3-11's V_G equation is an approximation (the slide says so), and is not used. |
| Brews, Solid-State Electronics 21, 345 (1978) | Primary paywalled. Galup-Montoro et al. 2007, eqs 5, 7, 9, 10a, 10b; the EE 316 slide 3-13 gives only the strong-inversion pinned form | PC, used as PD (16.18 (e)) | The charge sheet Q_I(φ_s) = −C_ox(V_G − V_FB − φ_s − γ√(φ_s − φ_t)) (eq 7), and the current as drift plus diffusion (eq 9) integrated between the source and drain surface potentials φ_s0 and φ_sL (eqs 10a, 10b). The surface potential comes from eq 5, which carries an author-dependent term χ(φ_s) "to circumvent some (numerical) problems … near the flat-band condition". |
| Brews, J. Appl. Phys. 45, 1276 (1974) (HF capacitance) | Abstract only (scilit) | PC on the premise; model PD | The HF capacitance of an inverted MOS-C keeps the total minority charge at its DC value, but lets the minority carriers move inside the inversion layer with a constant quasi-Fermi level across it. This is 16.2 (g)'s redistribution limit, which the plan derived independently. Brews reports a 7% error from neglecting that motion (strong inversion, 1e13 cm⁻³). Baccarani and Severi (IEEE, 1974) found ≤ 5% for practical dopings. |
| Sah, Noyce and Shockley, Proc. IRE 45, 1228 (1957) | Paywalled; abstract only | — (attribution only) | The reference does not use their closed-form approximations. It integrates U(ψ) of the spec's SRH over the QE profile, which is the definition (PD). |
| Shockley, BSTJ 28 (1949) | Read in TVR-0 | PA | Unchanged (I-9). |
| Wachutka (1990) | Not needed | — | TVB-28's reference is NiTCAD's own T2/T4 model equations (DECISIONS Unit 23), derived to S and Π here (PD), as DC7 did. 31.29 already carries the thermopower function (PD). |
| Nicollian and Brews, *MOS Physics and Technology*; Taur and Ning | Not reached | — | Not needed after Brews 1974's abstract and the two Pao–Sah reproductions; listed so they are not later cited as read. |

**Found in reading: the two Pao–Sah reproductions differ in the bulk term, and the difference is Galup-Montoro's
χ.**
- Reproduction (1) has the electron term (n_i²/N_a²)(e^{−V}(e^{ψ} − 1) − ψ), in units of V_T.
- Reproduction (2), eq 5 with χ = 0, has e^{−V}(e^{ψ} − ψ − 1).
- They differ by (n_i²/N_a²)(1 − e^{−V})ψ. For V > 0 that term makes E² negative for 0 < ψ < ψ* ≈ 2(n_i²/N_a²)(1 −
  e^{−V}) V_T, which is about 5e-16 V at N_a = 1e17. So form (1) is undefined near the bulk, and the inner integral
  of n/E would diverge logarithmically there.
- Physically, Pao–Sah's slice gives the electrons the channel's quasi-Fermi level V down to the bulk. In the bulk
  they are in fact at the substrate's level. Both the cutoff at ψ_B and χ are ways of not counting those electrons.
- NiTCAD solves the full problem, so the reference must bound the ambiguity, not choose between the forms. 16.18 (d)
  does that.

**16.18 Step 1: derivations (PD; each with the second method its code step will run at L0).**

*(a) TVB-13, the energy-balance boundary layer λ.* Linearize ED2, ED3 and Poisson about the uniform state (n = N,
field E₀, T_h with E_eff(T_h) = E₀, current J). Use modes e^{sx} with δJ = 0, unknown amplitudes (δψ, δn, δT), μ =
μ(T_h), μ′ = dμ/dT, κ = (5/2)(k²/q)Nμ T_h and G = (3/2)kN/τ_w:
- current: μqE₀ δn + μ′qE₀N δT − μNq s δψ + μk(T_h s δn + N s δT) = 0;
- Poisson: ε s² δψ = q δn;
- energy: (−κs² − (5/2)(k/q)J s + G) δT + J s δψ + (3/2)k(T_h − T_L)/τ_w δn = 0.

Eliminating δψ gives s · P(s) = 0, with the quartic
P(s) = [μkT_h s² + μqE₀ s − μNq²/ε](−κs² − (5/2)(k/q)Js + G) − (μ′qE₀N + μkNs)(Jq/ε + (3/2)k(T_h − T_L)s/τ_w).
- λ = 2/min|Re s| over P's four roots: the slowest mode, times 2 for safety.
- T_n is judged on the nodes with e^{−d/λ} ≤ 1e-13, d the distance to the nearer contact.
- The roots come from the quartic's companion matrix. The reference has no non-symmetric eigensolver, so the code
  step uses Durand–Kerner in dd, with its own L0 test.
- Second method: the linearized ODE is solved by Chebyshev collocation on a long interval, and its decay fitted.
- The signs follow ED2/ED3's conventions. The second method is what catches a sign slip, and it runs before the
  predicate is used.

*(b) TVB-10 HF, the redistribution limit.* Q_s(ψ_s, V) is TVB-09's KN charge with the electrons at a quasi-Fermi
level V that is uniform across the layer (E² as in reproduction (1), with V; holes at equilibrium). HF keeps
Q_n(ψ_s, V) at its DC value, so δV = −(∂Q_n/∂ψ_s / ∂Q_n/∂V) δψ_s, and
C_s,HF = −(∂Q_s/∂ψ_s + ∂Q_s/∂V · δV/δψ_s), C_HF = (1/C_ox + 1/C_s,HF)⁻¹.
- The partials come by complex step on the reference.
- Second method: a central difference in dd on Q_s along the constrained path, Q_n held fixed by Brent.
- Limits: in depletion C_s,HF → the depletion capacitance (Q_n → 0); with "V frozen" in place of "Q_n fixed" it is
  the simpler formula Brews improved on, and the reference's report prints the difference.
- ε_approx: the formula is exact for its limit. The departure from the limit is frequency-dependent, so it is measured
  as the spread of NiTCAD's C over the judged decade: ε_HF = max |C(f_i) − C(f_j)| over three frequencies spanning a
  decade, with the predicate ε_HF ≤ 1e-3 C (else NOT APPLICABLE, "not on the HF plateau").

*(c) TVB-15 A: the injection ladder (replaces a hand-set ε_approx; DT-Q14).* DD departs from the linear model by
O(Δn/N_A): the ambipolar field and the SRH nonlinearity, both first order. Rather than bound the constant, run DD at
G = g₀ × {1, 2, 4, 8} with g₀ giving max Δn/N_A = 1e-3. Then:
- the departure e(G) = ‖Δn_sim/G − Δn_lin/G‖∞ must have an observed order in G within [0.8, 1.3];
- its extrapolation to G → 0 must agree with the E comparison's limit within that comparison's U.

This measures the approximation error's order rather than assuming a constant. A departure of any other order is a
finding.

*(d) TVB-20 and 21: the channel-length limit (replaces a hand-set ε_approx; DT-Q13).* Pao–Sah is the L → ∞ limit of
the 2D problem at fixed W/L scaling: the gradual-channel error falls with L, and the source and drain ends add
ΔL/L. So run NiTCAD's long-channel MOSFET at L = L₀ × {1, 2, 4, 8} (L₀ ≥ 20 × the maximum depletion depth), with
Q(L) = I_D L/W. Then:
- the observed order of Q(L) in 1/L must lie in [0.8, 1.3], with the asymptotic check (order spread ≤ 0.1 on four
  lengths);
- Q_∞ by Richardson, with U_L = the GCI of the ladder;
- ε_approx = U_L + ε_cutoff, where ε_cutoff = |I_PS(cutoff ψ_B) − I_PS(cutoff ψ_B/2)| + |I_PS(form 1) − I_PS(form 2)|,
  from the reference itself (16.17's finding);
- applicability: V_G with ψ_s(V_S) ≥ ψ_B + 3V_T (weak inversion and above). Deeper subthreshold is NOT APPLICABLE,
  because there the bulk ambiguity is not small. Constant μ, no field mobility, no Lombardi, Boltzmann, no BGN.

S (TVB-21) is taken from the same ladder: d log₁₀ I/dV_G by a dd central difference on the reference, and the
ladder's Q_∞ at the V_G grid.
- Second method for the reference: Galup-Montoro's eq 4, g_d = −(μW/L)Q_I(V_D), against a dd difference of I_D in V_D.
  It is an identity of the double integral.

*(e) Brews' charge sheet (second A, diagnostic).* Integrating eq 9 with eq 7 between φ_s0 and φ_sL gives eqs 10a
and 10b. That integration is redone here (PD), and the second method is quadrature of eq 9 in dd. It is reported
beside Pao–Sah, as Brews himself compared it, and is never gated.

*(f) TVB-16 ε_approx (first-order terms, each measured or computed, summed as magnitudes per TVD-17).*
ε_J = J_ref × [ max(Δn/N) at the two depletion edges (measured on the solved state) + Δφ/V_T (the measured
quasi-Fermi drop across the depletion region) + 3L_D/(W − x_edge) (the edge's definition; short base only) ].
- Applicability: ε_J ≤ 0.1 J_ref (else NOT APPLICABLE).
- Metric: the log error, with q_floor the current resolution (Units 15 and 32).
- In (ii) the SNS integral is part of J_ref, and the reverse exponent is judged as d ln|J|/d ln(V_bi − V) against the
  reference's own exponent.

*(g) TVB-18 ε_approx.* ε_Y = |Y_ref| × [ max(Δn/N) + Δφ/V_T + J_SNS/J_diff (the depletion recombination eq 4.23
omits) + 3L_D/|L_eff(ω)| ], with L_eff = √(Dτ/(1 + iωτ)).
- Applicability: |L_eff(ω)| ≥ 30 L_D and ε_Y ≤ 0.1 |Y_ref|.
- Second method for the reference: the linearized minority ODE solved directly by complex Chebyshev collocation.

*(h) TVB-19 references (E).* Dielectric relaxation: Q(t) = Q₀e^{−t/τ_d} exactly for the linear majority equation.
Haynes–Shockley: Δp(x, t) = (G₀/√(4πDt)) exp(−(x − x₀ − μE t)²/(4Dt) − t/τ), convolved with the source's spatial
Gaussian (closed form: the variances add) and its time profile (Gauss–Kronrod in dd). Second method: the PDE by
method of lines with Dormand–Prince.

**16.19 Proposed A-class criteria (frozen on the owner's ruling, DT-Q13–Q17, before any code that runs them).**

| ID | Criterion | Applicability (else NOT APPLICABLE, predicate and value printed) |
|---|---|---|
| TVB-10 HF | \|C_sim − C_ref\| ≤ U_num + U_ref + ε_HF; ε_HF from 16.18 (b) | plateau spread ≤ 1e-3 C; strong inversion and depletion V_G grid |
| TVB-15 A | departure order in G ∈ [0.8, 1.3]; the G → 0 extrapolation agrees with the E comparison within its U | max Δn/N_A ≤ 1e-1 at the largest G |
| TVB-16 | log error ≤ log₁₀(1 + (U_num + U_ref + ε_J)/J_ref) | ε_J ≤ 0.1 J_ref; V ∈ [0.2, 0.6] V forward, [−5, −0.5] V reverse ((ii) only) |
| TVB-18 | \|Y_sim − Y_ref\| ≤ U_num + U_ref + ε_Y, magnitude and phase | 16.18 (g) |
| TVB-20 | the L-ladder order ∈ [0.8, 1.3], asymptotic check passed; \|Q_∞ − Q_PS\| ≤ U_num + U_ref + ε_approx; ε_approx ≤ 1e-2 Q_PS (else NOT EVALUATED) | 16.18 (d) |
| TVB-21 | \|S_∞ − S_PS\| ≤ U_num + U_ref + ε_approx (the same rule on the derivative) | as TVB-20, at V_D ≤ 3V_T |

**16.20 Owner questions from step 1 (recommendations in the last column).**

| # | Question | Recommendation |
|---|---|---|
| DT-Q13 | TVB-20/21: judge NiTCAD's channel-length limit (four lengths, Richardson) against Pao–Sah, with ε_approx = the ladder's GCI + the reference's own bulk-term ambiguity, instead of a hand-set gradual-channel bound? | Yes. It measures the approximation's order instead of assuming its size |
| DT-Q14 | TVB-15 A as an injection ladder (order in G, G → 0 limit = the E case) instead of a constant times Δn/N? | Yes, for the same reason |
| DT-Q15 | TVB-13: λ from the linearized quartic (16.18 (a)), × 2, with a Durand–Kerner root finder added to the reference (L0-tested) and the Chebyshev decay fit as the second method? | Yes |
| DT-Q16 | TVB-16/18: ε_approx as the sum of the named first-order terms of 16.18 (f), (g), each measured on the solved state or computed in the reference, with ε ≤ 0.1 of the reference as the applicability bound? | Yes. No source states a bound; each term is a named mechanism |
| DT-Q17 | TVB-10 HF: the redistribution limit (Brews 1974's premise, abstract only; derived here), ε from the measured plateau spread? | Yes. The "frozen V" variant printed as a diagnostic |
| DT-Q18 | Provenance: Pao–Sah PB (two reproductions agreeing on the double integral); Brews 1978 PC used only as PD with quadrature; SNS PD (the definition); Brews 1974 PD with the abstract recorded? | Yes |

Step 1 is complete apart from these rulings. Nothing in 16.8's frozen D/E criteria changed.

**16.21 Rulings on step 1 (owner "go" after the recommendations, 2026-11-06).** DT-Q13 through DT-Q18 are approved
as recommended. The A-class criteria of 16.19 are **frozen** from this commit, before any code that runs them (TVM2).

**16.22 Step 2, the reach of DT-Q2 A, measured before any change (16.2 (o)): it exceeds the ruled bound.** Removing
`device.cpp:75` and choosing an electrode node as the reference node is not enough. Four places in `assemble` take
a *semiconductor* from `Device::reference_node()`:
- `scaling.cpp:20`: the material behind the density and potential scales;
- `scaled_device.cpp:74`: the intrinsic-level depth that fixes the potential's reference;
- `gate_tunnels.cpp:64`: the same depth;
- `electrothermal.cpp:111`: κ_ref, by the reference node's region.

An insulator-only device needs a new source for each: for example a `DeviceDescription::reference_material`, a
semiconductor used only for scales, plus refusals of tunnelling and electrothermal. That touches `device` and three
or four `assemble` files. 16.2 (o) says to stop for the owner, so step 2 has stopped here. No code is written.

**A route that needs no change to NiTCAD (found while measuring).** `check_topology` (`device.cpp:281`) already
accepts a mesh of several connected parts: each part needs an ohmic contact or an electrode. So the benchmark device
can be the **capacitor component plus a disjoint silicon island** (a short n-type bar with one ohmic contact):
- The island supplies the reference node, and hence every scale, unchanged.
- No edge joins the two parts, so the box method's Poisson rows of the capacitor see nothing of the island. The
  electrodes' charges are exactly the capacitor's.
- Assembling the mesh: the two SimplexMeshes are concatenated with offset node ids (a test-side helper in `verify`,
  with its own test: each part's dual equals its own dual alone), then passed through `make_dual_mesh`.
  `make_dual_mesh` has no connectivity rule (`simplex_mesh.hpp:63`), and its P2 check only bites if the island sits
  inside a boundary facet's diametral circle, so the island is placed far away. `certify_interfaces` then certifies
  the union.
- In 1D, a graph mesh of the two chains with `certify_interfaces_1d`, which judges only edges joining two regions.
- TVB-06's frozen "2D tensor" family: each rectangle is split into two right triangles. Their circumcentres are the
  rectangle centres, so the dual is the tensor box method. The test asserts that equality before the benchmark uses it.
- Risk, to be measured in this route's first commit: whether the 2D generator accepts the coaxial annulus (sectors
  of one region whose union has a hole). If it does not, TVB-05 needs the owner again.

**DT-Q19.** (A) DT-Q2 A as ruled, widened: `reference_material` in `device` plus the `assemble` changes above, each
bit-identical for devices with a semiconductor, with their own tests. (A′) The disjoint-island fixture: no NiTCAD
change at all; the union helper and its test live in `verify`. (B) TVB-05/06 unbuildable. **Recommendation: A′.** It
is exact by construction (no shared edge), keeps NiTCAD untouched as CV-Q3 intended, and its one assumption (the parts
are independent) is tested directly rather than argued.

**DT-Q19 ruling: A** (owner, "A", 2026-11-06), not the recommended A′.

**16.23 Step 2 as built (DT-Q19 A).**

- `device`:
  - `DeviceDescription::reference_material` (optional). A description with no semiconductor region must name one, and its temperature is checked against it. A description with a semiconductor region must not name one.
  - `Device::has_semiconductor()` and `Device::reference_material()`. For a device with a semiconductor, `reference_material()` returns the same object as `material(reference_node())`, so every existing device is bit-identical by construction. A test asserts the address equality.
  - An insulator-only device's reference node is its lowest electrode node. Contact nodes are validated strictly increasing, so it is the first node of each electrode.
- `assemble`:
  - `scaling.cpp`, `scaled_device.cpp` and `gate_tunnels.cpp` read `reference_material()`.
  - `drift_diffusion` refuses the electrothermal model and any tunnelling contact on an insulator-only device.
- `io` (beyond DT-Q19's `device` + `assemble`; see DT-Q21): `save_device` refuses an insulator-only device, because the device-file format stores no reference material. Nothing is written.
- Tests:
  - `device/insulator_test`: the rules above.
  - `solve/insulator_test`: a 1D SiO₂ plate capacitor through `solve_small_signal`. Im Y/ω = ±C and Re Y ≈ 0 to 64 ε N², at 1 kHz and 1 GHz, at 0 V and 1 V. Both refusals are also checked.
  - `io/device_file_test`: the refusal.
- `verify/reference/capacitor`:
  - ε₀ in F/cm from CODATA 2018's printed digits (`codata2018::to_dd`, shift −2).
  - `series_capacitance` and `coaxial_capacitance` in double-double.
  - L0 tests: each form against a Gauss–Kronrod integral of its elastance density (≤ 16 and 64 u_dd), the one-layer and split-layer identities, and the thin-shell limit (≤ 1e-6).
- Mutation: 14 mutants (6 device, 2 assemble, 6 capacitor), all caught. The tunnelling-refusal mutant survived the first run, and the test was added before the rerun.
- Suites:
  - Release, the whole suite with the device change (`full` label and app tests excluded, −j 4): 171/171 pass.
  - Debug serial verify (`ctest -R verify -LE full`): 62/62 pass in 74.97 s and 76.82 s.
    - TVR-4's share is 0.13 s (`device_theory` 0.10 s, `capacitor` 0.03 s), within 16.7's 2 s.
    - The rest is the baseline: the unchanged tests ran 4–7 % above the 69.89 s of 15.20 (`fermi_dirac` +1.0 s, C-transient +1.1 s), as in step 8's 71–80 s under load.
  - The first Debug run had TVB-06's full case (2D 32²) under `[tvr4]`, costing 5 s. 16.7 puts it in `full`, so it moved to `[.tvr4_full]`.
- **TVB-06 (D): pass at every level**, against a criterion of 64 ε N² (1.5e-11 at N = 32, 2.3e-10 at N = 128):

  | Level | C (−Im Y_ab/ω) | C (Im Y_aa/ω) | E_ox/E_hk |
  |---|---|---|---|
  | 1D N = 32 | 2.9e-15 | 5.3e-15 | 7.9e-16 |
  | 1D N = 128 | 2.2e-14 | 2.4e-14 | 6.3e-16 |
  | 2D 32² | 3.6e-14 | 1.4e-15 | 7.9e-16 |

- **TVB-05 (E): the 16.22 risk did not occur.** `generate_mesh` accepts the annulus of M quadrilateral sectors.
  - Relative errors: 3.687e-5 at M = 256, 1.513e-4 at 128, 6.238e-4 at 64 and 2.729e-3 at 32. The backward errors are ≤ 7.2e-16.
  - (b) accuracy: **pass** (3.7e-5 ≤ 1e-3).
  - (a) order: the successive orders, finest pair first, are 2.0545, 2.0825 and 2.1846. p̂ = 2.054 is in the band [1.8, 2.3]. The spread is **0.1021**, above the asymptotic-range limit of 0.1 (§4.2, `metrics.hpp`), so `judge_order` returns **not_evaluated** ("not in the asymptotic range").
  - The excess comes from the coarsest pair. At M = 32, h = 1.3e-5 cm, which is a/8.
  - The frozen expected outcome (pass) is unmet, so `verify_benchmark.tvr4_full` fails (37 s in Release). Nothing has been changed after the run (TVM2).

**DT-Q20 (TVB-05's order).**
- (A) Record it as found: the order comparison is not_evaluated, with p̂ = 2.054, the spread of 0.102 and the coarsest pair named. The frozen expected outcome in the test becomes not_evaluated for that comparison, and the finest-accuracy pass stands. No criterion and no ladder changes.
- (B) Shift the ladder one level finer (M = 64…512), so that the four-level check excludes M = 32. This changes the frozen ladder (16.8) after a run, which TVM2 forbids in spirit. It also costs about 4× the finest level.
- (C) Keep the ladder and judge p̂ on the finest three levels (the 3-level rule of 14.7). This changes the method after the run, the same objection as B.

**Recommendation: A.** It is the only option that leaves 16.8 as frozen. The miss is 2 % of the limit on the coarsest pair, and the accuracy passes with a factor of 27 to spare.

**DT-Q21 (the `io` refusal, outside DT-Q19's `device` + `assemble`).**
- (A) Keep it. `save_device` refuses an insulator-only device, and the format is unchanged.
- (B) Widen the device-file format to store `reference_material`. This is a format change: the version, the reader, the writer and the round-trip tests.
- (C) Remove the refusal. `save_device` would then write a file whose `load_device` fails ("needs at least one semiconductor region or a reference_material").

**Recommendation: A.** It is the smallest change that keeps the save–load round trip honest. B can follow if a deck ever needs an insulator-only device.

**16.24 Rulings DT-Q20 and DT-Q21** (owner, "choose recommended", 2026-11-06).
- **DT-Q20: A.** TVB-05 is recorded as found.
  - The order comparison is not_evaluated: p̂ = 2.054 is in the band, but the spread is 0.1021, from the M = 32 pair.
  - The finest-accuracy comparison passes, and the overall outcome is not_evaluated.
  - The frozen expected outcome in `device_theory_test` now states this, so `verify_benchmark.tvr4_full` asserts it.
  - 16.8 (criteria and ladder) is unchanged.
- **DT-Q21: A.** `save_device` keeps refusing an insulator-only device, and the device-file format is unchanged.
  Widening the format (B) waits until a deck needs an insulator-only device.

**16.25 Step 3 as built: the hook's time profile (DT-Q3 A).**
- `assemble/models.hpp`:
  - `ManufacturedSources` gains `profile_times_s` (τ₀ < τ₁ < …) and `profile_weights` (one entry more than the times).
  - `manufactured_weight(ms, t)` gives w(t), piecewise constant and left-continuous: w₀ for t ≤ τ₀ (the steady state at t = 0 included), w_k on (τ_{k−1}, τ_k], and the last weight after the last time.
  - With both vectors empty, w ≡ 1.
- `DriftDiffusion::TimeStep` gains `time_s`, the step's end time, defaulting to 0. Existing aggregate initialisers are unchanged.
- The manufactured term multiplies each source by w at the step's end time, or by w(0) with no step. Without a profile w = 1, and 1·y = y exactly, so every result is bit-identical. A test asserts that a profile of ones is bit-identical to none.
- `create` refuses a malformed profile:
  - a weight count other than the time count + 1;
  - a time that is not finite, not positive, or not strictly increasing;
  - a weight that is not finite.
- `solve_transient` adds every τ_k in (0, t_end] to its restarts, so steps land on the jumps and restart with backward Euler, as at a waveform breakpoint. It passes t_new to each step. Each step therefore sees one weight, and a piecewise-constant source is integrated exactly.
- The profile is not in the deck vocabulary or the run record, as for TVQ3.
- Tests:
  - `assemble/manufactured_test`, three new cases: w(t) at and between the jumps; the term against w × the unweighted term at five times (1D and 2D), with no step meaning w(0); the bit-identity; and the refusals, 14 inputs.
  - `solve/transient_test`, one new case: a generation pulse G = 1e24 cm⁻³ s⁻¹ on (1, 2] ns in a 100 µm n-type bar (1e17, recombination off, fixed steps of 0.1 ns).
    - Steps land on both jumps and restart there.
    - Before the pulse the fields are the sourceless run's, bit for bit (10 snapshots).
    - At τ₁ the excess n and p in the middle are G·Δt = 1e15 to the printed digits (criterion 1e-3). At t_end they are still 1e15 (criterion 1e-2).
- Mutation: 16 mutants (2 models, 6 term, 5 validation, 3 transient), all caught. They include evaluating w at the step's start, dropping the restart, and `upper_bound` for `lower_bound`.
- Suite: Release, the whole suite (`full` label and app tests excluded, −j 4): 174/174 pass.

**16.26 Step 4 as built: minority_diffusion and TVB-15.**
- `verify/reference/minority_diffusion` (PD, base-free):
  - `srh_linear_inverse_lifetime`: 1/τ = (n₀ + p₀) / (τ_p(n₀ + n_ie) + τ_n(p₀ + n_ie)).
  - `minority_green`: the Dirichlet Green's function of −D u'' + u/τ, in an overflow-free exponential form.
  - `minority_excess`: the Green's-function integral of a Gaussian G, split at s = x, by qag15 in double-double to 0.5e-28.
  - `minority_collocation` (the second method): 8 pieces of degree 40, patched in u and u′, as a dense double-double LU.
- L0 tests (`verify.minority_diffusion`, 0.5 s in Debug):
  - the lifetime against a double-double difference quotient of the full SRH rate (≤ 1e-14);
  - the Green's function against its sinh form (≤ 64 u_dd), its symmetry and its derivative jump −1/D;
  - the integral against the infinite-domain erfc closed form on a bar with kL = 200 (≤ 1e-25);
  - the collocation against the integral (≤ 1e-22 of the peak), its tail estimate, and the throws.
- TVB-15 (full tier only, 16.7):
  - The bar: silicon, N_A = 1e16, L = 1 mm, ohmic ends, SRH only.
  - G is a Gaussian, σ = 100 µm, at L/2, through the TVQ3 hook.
  - Inputs (I1): D = μ_n V_T = 32.26 cm²/s and τ = 8.333e-6 s, giving a diffusion length of 164 µm.
  - The reference: I3 agreement 2.0e-28 of the peak; U_ref 2.0e-28.
  - **E (16.8): pass.** p̂ = 1.9993 (orders 1.9993, 1.9987, 1.9952; spread 0.0035). The finest error (N = 512) is 8.14e-6 against 1e-4. Max Δn/N_A = 5.0e-9 on every level (predicate 1e-8). U_iter ≤ 9.5e-16.
  - **A (16.19): pass.** The injection ladder at N = 512, g₀ for a peak of 1e-3 N_A:
    - the departure's order in G is 0.995 (then 0.990), in [0.8, 1.3];
    - the G → 0 limit by Richardson differs from the E solve by 2.93e-6, within U_E = 1.02e-5 (U_disc by the finest pair's GCI, F_s = 1.25, + U_iter + U_ref);
    - max Δn/N_A at 8 g₀ is 8.0e-3 (applicability 1e-1);
    - the ladder's U_iter is ≤ 2.9e-15, below a tenth of the smallest difference (CV-Q8, applied to the ladder).
  - The A record states its approximation error as measured (the `contract` rule for A).
- Construction change before commit (no criterion or ladder touched, TVM2):
  - The first build used L = 100 µm and σ = 10 µm. It passed all four comparisons (E 6.0e-6, p̂ 1.999; A order 0.997, limit 2.8e-6 within 7.6e-6).
  - A mutation then showed that a bar shorter than the diffusion length barely depends on τ.
  - The bar was lengthened to the spec's "long" (about 6 diffusion lengths, the same σ/L) and the run above was repeated.
- Mutation:
  - reference: 9 mutants, all caught. The tail-estimate mutant survived until the test asserted the estimate is computed.
  - benchmark: 13 mutants, all caught, including τ × 1.01, D × 1.01, G × 1.001, the predicates, the limit bound, the ladder's CV-Q8 rule and the U_ref gate.
  - One mutant, τ = τ_n in place of the linearization, survives. It is a relative 1e-6 change (n_ie/p₀), below the 1e-4 tolerance, so the benchmark cannot resolve it. This is recorded, not hidden.
- Timing and suites:
  - Release `verify_benchmark.tvr4_full` (TVB-05, TVB-06 2D, TVB-15) passes in 35.6 s, inside 16.7's 120 s.
  - Debug serial verify passes, 63/63 in 71.7 s.
  - `registry`: 44 benchmarks, TVB-15 in `full`.

**16.27 Step 5 as built: qe_junction and TVB-17.**
- `verify/reference/qe_junction` (PD, base-free):
  - The QE junction at bias V is TVR-2's equilibrium junction with n_i → m = n_i e^{V/2V_T}, the potential shifted by V/2. Bias enters only through m. Derived in the header; the L0 test checks it against `pn_junction` at m.
  - Q = ε E_max = √(2 q ε V_T Φ_n(u_j)), with u_j by Brent in double-double.
  - C = −dQ/dV by one dual-number Newton step from the converged root. It is exact to rounding (du_j/dv = −F_v/F_u).
  - The second method is complex step (h = 1e-200) on Q(v) with complex Newton for u_j.
  - **Plan deviation (method only):** 16.3 says "C by complex step on Q(V)", which is the second method here. Complex step cannot run through a double-double Brent root, so the primary is the dual number.
  - A guard refuses m < 1e-140 cm⁻³ (about −15 V at 300 K), where m² leaves the normal range.
- Derived (PD, recorded because 16.3 expected a reportable "ε/W difference"):
  - Dropping the junction's minority terms and taking n_n = N_D, p_p = N_A turns the first integral into exactly ε/W with V_bi − V − 2V_T.
  - So ε/W(−2V_T) equals the QE C up to terms exponentially small in the barrier and m²/N².
  - Measured: 9.1e-8 at −0.5 V, 1.4e-10 at −1 V, 4.6e-16 at −2 V.
- L0 tests (`verify.qe_junction`):
  - zero bias against `pn_junction` (≤ 1e-28);
  - the shift identity at four biases (≤ 1e-15, the double rounding of m);
  - dual C against complex step (≤ 1e-13) and against a Richardson central difference with exact 2⁻²⁰ steps (≤ 1e-14);
  - ε/W converging to C;
  - the range guard and the input errors.
- Mutation (reference): 11 mutants, 10 caught. The survivor is an equivalent mutant: the derivative part of F_u, which drops out of the Newton step's derivative at the root. The range-guard mutant survived until a −18 V case refused only by the guard was added.
- TVB-17 (full tier only):
  - The diode: TVB-07's junction on [−1.5, 1.5] µm, no G−R, the anode on the p side.
  - V = −0.5, −1, −2 V; C = Im Y_aa/ω at 1 kHz; N = 64…512.
  - DT-Q4: applicable at every level, with quasi-Fermi drops ≤ 5.6e-11 V against 1e-3 V_T = 2.6e-5 V.
  - U_ref ≤ 1e-11: dual vs complex step ≤ 4.6e-16, plus the contact tail of 1.9e-13, 8.6e-13 and 9.9e-12.
  - U_iter ≤ 5.9e-9, below U_disc/10 everywhere.

  | V | finest error (N = 512) | (b) ≤ 1e-3 | orders, finest pair first | spread | (a) |
  |---|---|---|---|---|---|
  | −0.5 V | 5.65e-5 | pass | 1.9986, 2.0010, 2.2256 | 0.225 | not_evaluated |
  | −1 V | 3.96e-5 | pass | 2.0003, 2.0011, 1.8184 | 0.183 | not_evaluated |
  | −2 V | 2.48e-5 | pass | 2.0002, 2.0007, 1.8675 | 0.133 | not_evaluated |

  - The coarsest pair (N = 64, h = 47.6 nm, about the 41 nm p-side Debye length) puts the spread above the asymptotic limit of 0.1. The finer pairs agree with 2 to 1e-3.
  - The frozen expected outcome (pass) is unmet on the three order comparisons, so `verify_benchmark.tvr4_full` fails (36.6 s). Nothing has been changed after the run (TVM2).
  - The benchmark's own mutation run waits for the ruling: with the test failing, every mutant would look caught.
- Timing: Debug serial verify passes, 64/64 in 72.5 s (`qe_junction` 0.02 s). `registry`: 45 benchmarks, TVB-17 in `full`.

**DT-Q22 (TVB-17's order).**
- (A) Record it as found, as DT-Q20 A: the three order comparisons are not_evaluated, with the orders and spreads above, and the three finest-accuracy passes stand. No criterion, ladder or construction changes.
- (B) Shorten the diode to TVB-07's L = 1 µm per side, so that h(64) = 32 nm. This changes my construction after a run, to make a check pass, so it is TVM2's concern again. The contact tail would rise to about 3e-6 at −2 V, still within U_ref ≤ 1e-4. Whether the spread would then meet 0.1 is not known.
- (C) Judge p̂ on the finest three levels (14.7's 3-level rule). This changes the method after the run.

**Recommendation: A.** It matches DT-Q20 and leaves 16.8 untouched. The finest pairs show the scheme's order cleanly at every bias.

**16.28 Ruling DT-Q22** (owner, "choose recommended", 2026-11-06).
- **DT-Q22: A.** TVB-17 is recorded as found, as DT-Q20 A.
  - At each of V = −0.5, −1 and −2 V, the C_j order comparison is not_evaluated ("not in the asymptotic range"): p̂ is 1.9986, 2.0003 and 2.0002 (finest pair), in the band, with spreads 0.225, 0.183 and 0.133 from the N = 64 pair.
  - The three finest-accuracy comparisons pass (5.65e-5, 3.96e-5, 2.48e-5 against ≤ 1e-3), and the overall outcome is not_evaluated.
  - The frozen expected outcome in `device_theory_test` now states this (each order not_evaluated with p̂ in [1.8, 2.3], each finest pass), so `verify_benchmark.tvr4_full` asserts it and passes.
  - 16.8 (criteria and ladder) and TVB-17's construction (16.27) are unchanged.
- Mutation (TVB-17 benchmark, against the new expected outcome): 12 of 12 caught. **Invalid, corrected in 16.30:** the filter matched no test; the rerun gives 11 of 12. The mutants were the DT-Q4 thresholds (electron and hole), the density cut, the U_ref gate, C_ref's sign, the tolerance, a 1 % scale on C, h, the frequency (1e11 Hz), N_D, p_th and the finest-accuracy judgement.
- `verify_benchmark.tvr4_full` passes.

**16.29 Found before step 6's code: two gaps (owner questions DT-Q23, DT-Q24).**
- **(a) The (h, tol) factorial was left out of steps 2, 4 and 5.**
  - 16.2 (p), approved with the plan (16.15, 16.16), requires every E and A benchmark to run the 2×2 (h, tol) factorial at its finest level and report ε_int. An ε_int above 10 % of U blocks an E pass (§4.1).
  - TVB-05, TVB-15 (E and A) and TVB-17 do not run it. `reference::interaction_contrast` and `combine_uncertainty` exist (metrics.hpp), but no benchmark calls them. Their recorded passes (16.23, 16.26, 16.28) therefore lack a term the plan requires.
  - This is my omission, not a finding about NiTCAD. Nothing has been changed.
- **(b) TVB-16's ideality and reverse exponent have no frozen tolerance.** 16.19 freezes only the log-error rule on J(V). 16.18 (f) says the reverse exponent "is judged" against the reference's own exponent, and 16.4 lists ideality, but neither has a criterion.

**DT-Q23 (the factorial).**
- (A) Before step 6, add it to TVB-05, TVB-15 and TVB-17 in one fix commit, and from now on to every E/A benchmark.
  - At the finest level: e(h, tol), e(h, tol/100), e(h/2, tol) and e(h/2, tol/100), with tol the Newton `tol_update` (1e-12 → 1e-14) and h/2 one more level.
  - ε_int is reported, and enters U by §4.1. If it exceeds 10 % of U, the E pass is blocked.
  - The frozen expected outcomes would change only if a pass is blocked, and that would be recorded as found.
  - If tol/100 cannot be reached (Newton stalls at rounding), that level is recorded as not_evaluated with the measured update, not loosened.
- (B) Defer it to the close-out (step 14), for all TVR-4 benchmarks at once.

**Recommendation: A.** The rule applies to each benchmark. Building TVB-16 onward with it, while 05/15/17 lack it, would leave the unit inconsistent.

**DT-Q24 (TVB-16's derivative quantities).**
- (A) Judge the ideality n = (1/V_T) dV/d ln J (forward) and the reverse exponent d ln|J|/d ln(V_bi − V) ((ii) only) by 16.19's rule applied to the derivative, as TVB-21 does for S:
  - |x_sim − x_ref| ≤ U_num + U_ref + ε_x;
  - U_num is the GCI of x on the ladder, and U_ref the reference's second-method difference;
  - ε_x carries ε_J from the two biases of the difference: ε_x = |x_ref| (ε_J,k/J_k + ε_J,k+1/J_k+1) / |Δ ln y|, with y = J for n and y = V_bi − V for the exponent;
  - applicability is as for J at both biases.
- (B) Report n and the exponent as diagnostics only. J(V) at each bias is already gated by the log error.

**Recommendation: A.** 16.18 (f) says the exponent is judged, and this rule is the one already frozen for S. ε_x is conservative (a sum of magnitudes).

**16.30 Rulings DT-Q23 and DT-Q24** (owner, "choose recommended for commercial TCAD", 2026-11-06), and the factorial fix as built.
- **DT-Q23: A.** The (h, tol) factorial is added to TVB-05, TVB-15 and TVB-17 before step 6, and is required of every later E/A benchmark.
- **DT-Q24: A.** TVB-16's ideality and reverse exponent are judged by 16.19's rule applied to the derivative (16.29).
- **The helper** (`benchmark::factorial`, `block_if_interacting`, registry.hpp) wraps `reference::interaction_contrast` and `combine_uncertainty`.
  - An E pass becomes not_evaluated with interacting errors, or when the factorial was not evaluated (4.3).
  - Fail, finding and not_applicable outcomes are kept.
  - L0 tests are in `registry_test`.
- **Results.** tol is the Newton `tol_update`, 1e-12 → 1e-14. U at (h, tol) = U_disc (GCI of the finest pair, F_s 1.25) + U_iter + U_ref.

  | Benchmark | levels (h, h/2) | ε_int | U | interacting |
  |---|---|---|---|---|
  | TVB-17, −0.5 V | N = 512, 1024 | 8.0e-13 | 7.06e-5 | no |
  | TVB-17, −1 V | N = 512, 1024 | −1.3e-12 | 4.95e-5 | no |
  | TVB-17, −2 V | N = 512, 1024 | −1.7e-12 | 3.11e-5 | no |
  | TVB-15 E | N = 512, 1023 | 0 (tol/100 rounds to the same state) | 1.02e-5 | no |
  | TVB-05 | M = 256, 512 | 0 (the AC solve is linear and direct) | 4.61e-5 | no |

  - TVB-15's A limit is now judged against U_E with |ε_int| included. Its frozen outcomes are unchanged.
  - TVB-17's frozen outcome (16.28) is unchanged.
- **TVB-05's h/2 level is beyond 16.7's CI budget.**
  - At M = 512 (131k nodes), mesh generation takes 70 s and each solve about 190 s. tvr4_full is 45 s locally without it.
  - Under 16.7, it runs in the `nightly` tier: the local `verify_full.tvb05`, registered with NITCAD_VERIFY_FULL_TIER. There the frozen outcome is: order not_evaluated (16.24), finest pass, ε_int reported and not interacting.
  - In `full` (and so in tvr4_full on CI), the factorial is not evaluated. TVB-05's finest accuracy is therefore held at not_evaluated: "the (h, tol) factorial was not evaluated". Its error, 3.7e-5, is still printed and asserted ≤ 1e-3.
  - This is the one frozen outcome that changes: TVB-05's finest pass in `full` becomes not_evaluated. TVB-05's overall outcome was already not_evaluated.
- **Found in my own mutation runs: 16.28's "12 of 12 caught" was invalid.**
  - Those mutants were run with the test's full name as the Catch2 filter. That name contains a comma, which Catch2 reads as a separator, so no test matched. The binary then exited non-zero, and every mutant counted as caught.
  - Rerun with a comma-free filter (`TVB-17*`): 11 of 12 benchmark mutants are caught. The survivor is p_th = 1: the order is never judged (16.28), so p_th cannot be observed.
  - Of the 4 factorial-wiring mutants, 3 are caught. The survivor inflates U_disc 1e9-fold: with no interaction, U's size cannot be observed.
  - No other mutation run in TVR-4 used a filter with a comma. Those runs used tags or comma-free wildcards, checked against this session's transcript.
  - TVB-15's 3 wiring mutants are all caught.
  - In the helper, 6 of 7 mutants are caught after one added check (a not_applicable outcome kept). The survivor, `ok = value == value`, is equivalent, since the value is finite whenever the status is ok.
- **Checks.**
  - `verify_benchmark.tvr4_full` passes in 45.2 s (36.6 s before), and `verify_full.tvb05` in 652 s (local; run alongside other work).
  - `[registry]` passes. Debug serial verify passes 64/64 in 75.5, 74.3 and 74.3 s (72.5 s at step 5). The fast tier is unchanged, and the slowest tests (fermi_dirac 21.6 s, consistency.transient 15.3 s) are untouched. So this is the machine's spread, but the 75 s budget now has under 1 s of margin, for the close-out (step 14).

**16.31 Step 6a as built: the diode reference (`verify/reference/diode`); DT-Q25 open before TVB-16's code.**
- **`diode_current(in, V)`** builds on TVR-2's equilibrium junction at m = n_i e^(v/2) (qe_junction's shift identity):
  - edges x_p = Q/(qN_A) and x_n = Q/(qN_D) from pn_junction at m, and bases W′ = L − x;
  - Shockley's finite-base law q(e^v − 1)[D n₀/L coth(W′/L)] per side, with L = √(Dτ) and 1/τ SRH's exact linearization at the bulk. Without SRH it is the short-base law q(e^v − 1) D n₀/W′;
  - the SNS current q∫U dx over [−x_p, x_n] on the reference's QE profile, as ∫U(u)/|u′(u)| du by Gauss–Kronrod in u, each side with its own lifetimes;
  - second method: tanh-sinh on the same integrals (16.3). The change of variable is checked at L0 against U(x) integrated in x with u(x) from pn_potential.
- **L0 (`diode_test`)**:
  - no SRH: the short-base law with edges recomputed from qe_junction's charge, to 1e-13;
  - long base (50 L) and short base (W′/L ~ 1e-6) limits of the coth law;
  - L against √(Dτ_n) on the p side and √(Dτ_p) on the n side;
  - the reported bulk minorities; SNS: the two quadratures agree; the sign; the exact 1/τ scaling;
  - at −5 V: below the bound q n_i(x_p + x_n)/(τ_n + τ_p), and within 10 % (measured 3.9 %) of the depletion-approximation estimate q n_i ΣW_eff/(τ_n + τ_p);
  - the x-space integral agrees to 1e-8;
  - zero current at zero bias, and the refusals.
- **Found while building (my errors, fixed before commit):**
  - (a) The n side's SRH linearization got its electron and hole densities swapped, so the hole lifetime came out as about τ_n. The added n-side length check catches it.
  - (b) A first second method integrated in x with a pn_potential solve per abscissa. It took 40 s in Release and 275 s in Debug, and passed vacuously: tanh-sinh's error estimate is absolute and capped at 1, meant for integrals of order one (tanh_sinh.hpp). Unscaled integrals of about 1e10 stop it at level 3, which is what happened.
    - Now every tanh-sinh integrand is scaled to order one.
    - The x-space check (2.2 s in Release, 15 s in Debug) is hidden ([.diode_profile]) and run by the Release-only ctest `verify.diode_profile`, to keep within CV-Q19's Debug budget.
    - The SNS L0 case runs two biases (0.86 s → about 0.4 s in Debug).
- Mutation: 20 of 20 caught. One survivor (the reported n_p0) was caught after a check was added.
- Release verify passes 68/68 (`verify.diode_profile` 2.1 s). Debug serial verify passes 65/65 in 73.0 s (`verify.diode` 0.70 s).

**TVB-16's construction, planned (my choices; nothing has run).**
- (i) Silicon, N_A = 1e16 for x < 0 and N_D = 2e16 for x > 0, 5 µm per side, so that 3L_D/W′ ≈ 0.025.
  - SRH, Auger and radiative off; Boltzmann statistics; no band-gap narrowing; no field mobility.
  - Caughey–Thomas mobility for both carriers. Holes are set to it explicitly, since the default is automatic Masetti.
  - Uniform nodes N = 128, 256, 512, 1024 (N even, so x = 0 is a face), with the factorial at N = 1024 and 2048.
  - Forward bias 0.2–0.6 V in steps of 0.1 V.
- (ii) short: the same device with SRH on (Scharfetter lifetimes at the local doping). Forward as (i), and reverse at −0.5, −1, −2 and −5 V.
- (ii) long: 500 µm per side, about 3 L_n and 10 L_p.
  - Nodes x = L sinh(βξ)/sinh β with ξ uniform and β = 7.5: the junction spacing is 8 nm at N = 1024. Same ladder.
- J_sim is the weighted terminal current (Unit 32), and its resolution goes into U_num.
  - U_num = GCI of J (Roache's three-grid form with the measured ratios, since h(N) ∝ 1/(N − 1) is not exactly constant-ratio) + the resolution + |ε_int|.
  - U_ref = the SNS second-method difference + the quadrature estimate.
- The inputs are NiTCAD's material data (I1), printed as the spec values (R-7): D = μ V_T with Caughey–Thomas μ at the side's doping, and Scharfetter τ.
- The legacy short-base value 1.321e-2 (DT-Q5) is printed as a cross-check, not gated.

**DT-Q25 (16.18 (f)'s quasi-Fermi term in reverse bias).**
- 16.18 (f) writes the term as Δφ/V_T, "the measured quasi-Fermi drop across the depletion region", as a first-order term.
  - In forward bias, J ∝ e^{(φ_p − φ_n)/V_T} at the edges, so the term is Δφ/V_T.
  - In reverse, the minority carriers' quasi-Fermi levels fall by about |V| across the region where their density is negligible. Read literally, the term is then about 20–200, so every reverse point would be NOT APPLICABLE and the reverse exponent, which 16.18 (f) says is judged, could not be.
  - The first-order sensitivity of both J_diff and J_SNS to a drop Δφ is e^v/|e^v − 1| × Δφ/V_T. This holds because both are ∝ (np − n_i²) at the edge or in the region, and np carries the e^{Δφ/V_T}.
- (A) Use ε_φ = e^v/|e^v − 1| × Δφ/V_T.
  - Δφ is the larger of |φ_n(x_n) − φ_n(−x_p)| and |φ_p(x_n) − φ_p(−x_p)|, interpolated at the reference's edges.
  - In forward bias this is Δφ/V_T to within e^{−v} (4.6e-4 at 0.2 V). In reverse it is ≈ e^v Δφ/V_T.
- (B) The literal Δφ/V_T: all reverse points NOT APPLICABLE, and the reverse exponent recorded as not judgeable.
- (C) Measure the drop only where the carrier's density exceeds 1e-6 of the local doping (DT-Q4's predicate). In forward bias this would exclude the minority edge itself (n(−x_p) ≈ 2e7 cm⁻³ at 0.2 V, against 1e10) and understate the term where it was written.

**Recommendation: A.** It is the first-order term 16.18 (f) asks for, and equals the frozen text in forward bias, where that text applies. Nothing else in 16.19 changes.

**16.32 Ruling DT-Q25** (owner, "choose recommended and write", 2026-11-06), **and step 6b as built: TVB-16. DT-Q26 and DT-Q27 open.**
- **DT-Q25: A.** ε_φ = e^v/|e^v − 1| × Δφ/V_T, with Δφ the larger quasi-Fermi drop between the reference's edges.
- **As built**, as planned in 16.31:
  - `tvb16.cpp`: three configurations, the ladder N = 128 … 1024, and the factorial at N = 1024/2048 with tol 1e-12/1e-14.
  - DT-Q24's derivative comparisons: 4 idealities per configuration, and 3 reverse exponents per (ii) configuration. 41 comparisons in all.
  - Registered in `full` (46 benchmarks). The case is in `device_theory_test` [.tvr4_full], with the frozen expected outcome "every comparison passes". It runs in 9.7 s (Release).
  - **Found while running (my construction, not NiTCAD).** A factorial tol/100 state that Newton cannot reach no longer aborts the run. Per DT-Q23 A it is recorded with its last update, and the comparison is not_evaluated. Newton stalls at 1e-14 at −5 V (last updates 1.4e-14 and 7.3e-14 to 2.0e-13) and at 0.5 V in (ii) short.
- **Results** (first run, nothing changed after it): 20 pass, 15 not_evaluated, 6 fail. The frozen outcome is unmet, so `verify_benchmark.tvr4_full` fails.
  - (i) no G−R, short base: all 9 pass. Log errors are 2.4e-6 to 6.4e-3, against bounds of 1.3e-2 to 2.1e-2. The ideality errors are 2.5e-6 to 3.7e-3.
  - (ii) short:
    - J passes at 0.4–0.6 V and at −0.5 to −2 V.
    - not_evaluated at 0.2 and 0.3 V: the three-grid differences oscillate, so there is no GCI.
    - not_evaluated at −5 V: the factorial's tol/100 was not reached.
    - The derivatives pass where both of their J values are estimated.
  - (ii) long: the 6 fails are J at 0.5 and 0.6 V and at −1 and −2 V, the 0.5–0.6 V ideality, and the −1 to −2 V exponent. They come from two causes:
    - **(a) Ohmic drop in the 500 µm neutral bases (forward bias).**
      - The deviation is converged in the mesh: e(h) and e(h/2) are −2.945e-3 and −2.928e-3 at 0.5 V, and −12.383 % and −12.384 % at 0.6 V. It is not discretization.
      - The majority current crosses both bases, whose series resistance is R = Σ W′/(q N μ_maj) ≈ 0.05 × (1.54 + 0.27) ≈ 0.090 Ω·cm². At 0.6 V the estimate J R / V_T ≈ 0.17 gives about −15 %, against −12.4 % measured; at 0.5 V it gives 0.40 %, against 0.29 %. Near the junction part of the current is minority current, so the true drop is smaller than this estimate.
      - Shockley's law assumes the whole of V falls across the junction. 16.18 (f)'s three ε_J terms have no series-resistance term. In the 5 µm devices the term is about 2e-4 of J at 0.6 V, which is why they pass.
    - **(b) U_disc is underestimated (reverse bias).**
      - At −1 and −2 V the error falls as h²: e(h/2)/e(h) = 0.249 in the factorial, so J_sim tends to J_ref.
      - But the three-grid order on N = 1024/512/256 of the graded mesh is 3.80 and 2.71, and the GCI with that order is about 4× smaller than the error it should bound.
      - §4.2 says an order outside the asymptotic range "is not reported as an order". The benchmark computes the GCI without that check.
- The DT-Q5 cross-check, printed only: the W(V) short-base law on bias:144's inputs gives 1.2786e-2 A/cm² at 0.5 V, against its W(0) value of 1.321e-2. bias:144's legacy pin is 1.280e-2.
- The TVB-16 benchmark mutation run waits for the rulings, as for DT-Q22.

**DT-Q26 (U_disc outside the asymptotic range).**
- (A) Apply §4.2's asymptotic check (successive orders within 0.1 on the 4 levels) to every J and derivative ladder.
  - Outside the range, the order is not reported, U_num is "not estimated" and the comparison is not_evaluated (§4.3).
  - This is §4.2's existing rule, which the code left out. How many comparisons it leaves judged has not been measured.
- (B) Outside the range, use Roache's two-grid GCI with the assumed order: F_s = 3.0 and p = p_th = 2 on the finest pair, N = 2048 and 1024 (already solved for the factorial).
  - This changes 16.31's U_num definition after a run.
- (C) Keep it as run: the −1 and −2 V long-base fails stand, as found.

**Recommendation: A.** It is the rule §4.2 already states, and it introduces no new choice.

**DT-Q27 (the series resistance missing from ε_J).**
- (A) Record it as found. 16.18 (f)'s ε_J omits the neutral bases' ohmic drop, which dominates the long base at ≥ 0.5 V.
  - The three forward comparisons it affects (J at 0.5 and 0.6 V, the 0.5–0.6 V ideality) become findings, through the known-findings table citing this entry, with the estimate above printed.
  - This is about the reference's approximation, not NiTCAD.
- (B) Add a fourth named first-order term, ε_R = J_ref R_base / V_T, with R_base computed from the spec, and rerun.
  - This adds to a frozen criterion after a run (TVM2), so it needs your explicit authorization.
- (C) Make the long base's forward points NOT APPLICABLE where J R_base / V_T > 1e-3 (a new predicate). Also a change after the run.

**Recommendation: A** for TVR-4: the frozen criterion stands, and the omission is recorded with its measured size. B would be the right design for a later revision of 16.18 (f), as its own ruling before any run.

**16.33 Rulings DT-Q26, DT-Q27** (owner, "choose recommended", 2026-11-06), **as applied; TVB-16 rerun. DT-Q28 open.**
- **DT-Q26: A.** `tvb16.cpp` now has `disc()`. It computes the three-grid GCI on levels 1024/512/256 (order p012) and on 512/256/128 (order p123). U_disc is estimated only when both orders exist and |p012 − p123| ≤ 0.1 (`reference::asymptotic_order_spread`, §4.2). Otherwise the comparison is not_evaluated (§4.3), and its reason says which case applies: oscillatory or non-positive finest order, no order on the coarser triple, or "outside the asymptotic range (orders a, b)". This applies to the J ladders and to the derivative ladders. Both orders are printed.
- **DT-Q27: A.** `registry.cpp` gets known finding **(x)**: TVB-16, quantity prefixes "J (ii) SRH, long base" and "ideality (ii) SRH, long base", region V ≥ 0.5 V (for the ideality, its lower bias). The ruling text cites this entry. The labels (viii) and (ix) were already used by TVB-31.12 and 31.15 (MI-Q15), so this one is (x). There are now 27 entries.
- **Found while applying (my construction, from 8ce7fe9).** `registry_test` was not updated when TVB-16 was registered. Its id list had 45 names for 46 benchmarks: from TVB-16 on, the checks compared the wrong names, `ids[45]` read past the array, and TVB-16 was not in the full-tier set. So 35 checks failed in `verify_benchmark`'s default run. 16.32's "Release verify 68/68" was therefore wrong for that commit. Fixed:
  - TVB-16 is added to the list and to the full-tier set;
  - a new check `REQUIRE(std::size(ids) == b.size())`.
- **Results** (Release, first run after the rulings, nothing changed after it): **1 pass, 39 not_evaluated, 1 fail**. `tvr4_full` still fails.
  - **pass:** J of (i) at 0.2 V only. Its orders are 1.949 and 1.983, and every other (i) point is outside the range. For example, J at 0.5 V has orders 1.81 and 0.34, and the 0.2–0.3 V ideality has 2.05 and 2.92.
  - **(ii) short base, all 17 not_evaluated.** At 0.2 and 0.3 V the finest triple oscillates. At 0.4–0.6 V and −1 and −5 V the coarser triple has no order. At −0.5 and −2 V the orders differ by more than 0.1 (4.09/3.32 and 2.69/3.84).
  - **(ii) long base.** Finding (x) changes nothing in this run: J at 0.5 and 0.6 V and the 0.5–0.6 V ideality are now not_evaluated (orders 2.09/7.99, 1.06/0.35, 0.77/6.38), so no fail remains for it to classify. The one **fail** is J at −2 V:
    - orders 2.713 and 2.797, inside 0.1, so U_disc is estimated: U_num/J = 1.15e-4;
    - but the error is 1.705e-4 (log error 7.41e-5 against a bound of 5.01e-5);
    - the factorial shows e(h/2)/e(h) = 0.249, so between N = 1024 and 2048 the error falls as h² and is discretization error;
    - the two three-grid orders agree with each other, but both are pre-asymptotic, near 2.7 against the measured 2.0. §4.2's check does not detect this, and the GCI is about 1.5 times too small.
  - **The derivatives:** every one is not_evaluated. Each needs its own ladder in range and both of its J values estimated.
- **Diagnostic** (a temporary build, reverted, not committed): with the factorial's N = 2048 state as the finest level, the J ladders 2048/1024/512/256/128 are within 0.1 at only 4 of 23 points: (i) at 0.2, 0.3 and 0.4 V, and the long base at −5 V. So a one-level-finer ladder does not bring the graded mesh into the asymptotic range.
- The TVB-16 mutation run and the frozen expected outcome wait for DT-Q28.

**DT-Q28 (TVB-16's ladder is almost never in the asymptotic range).** Under DT-Q26 A, 39 of 41 comparisons are not_evaluated. The one judged fail is a GCI that is too small, at orders that pass the check.
- (A) **Record it as found.**
  - The frozen expected outcome becomes the found outcomes, listed per comparison.
  - The −2 V long-base J becomes a known finding (xi): "the three-grid GCI underestimates on the graded mesh at in-range but pre-asymptotic orders (2.71/2.80 against the measured 2.0)".
  - TVB-16 then judges one point. Its test passes, but in practice it verifies almost nothing.
- (B) **Roache's two-grid estimate as the fallback.**
  - Where the three-grid estimate is unavailable or out of range, U_disc of the N = 1024 value is F_s |J_2048 − J_1024| r^p/(r^p − 1), with F_s = 3 and p = p_th = 2 (the coarse-solution form; the N = 2048 states are already solved for the factorial). Where the three-grid estimate is in range, it is kept, and the larger of the two is used.
  - It is DT-Q26's option B, now with A's result measured. It widens U_num after a run (TVM2), so it needs your explicit authorization.
  - Effect, from numbers already printed: at −2 V long, 3 × 1.28e-4 × 4/3 = 5.1e-4 ≥ 1.705e-4, so the fail would pass. Its effect elsewhere has not been computed.
  - p_th = 2 is supported where the factorial shows it (0.249 at −1 and −2 V). Where e(h/2)/e(h) is outside 0.25 ± 0.05 (a band this ruling would fix before any run), the fallback is not estimated either.
- (C) **A finer, longer ladder** (N = 4096 … 512) and rerun. This changes 16.31's frozen ladder after a run, costs about 4× the 9.7 s, and the diagnostic above gives little reason to expect the range.

**Recommendation: B.** A benchmark that judges 1 of 41 points does not verify the diode. Roache's F_s = 3 two-grid estimate with an assumed order is the documented practice when the observed order cannot be trusted (Roache, Verification and Validation in Computational Science and Engineering, 1998: F_s = 3 for two grids). The condition that e(h/2)/e(h) be near 1/4 keeps it from being applied blindly. It is a post-run change, so it needs your "authorize". A is the strict reading of TVM2.

**16.34 Step 7a as built: the admittance reference (`verify/reference/admittance`); DT-Q29–DT-Q31 open before TVB-18's code.** TVB-18's code also waits for DT-Q28, because its U_num follows the same GCI rule.
- **`diode_admittance(in, V, ω list)`** builds on diode.hpp's operating point (the QE edges at V, the bases W′, U₀ = n₀(e^v − 1)). Per base, with κ = √((1/τ + iω)/D) (1/τ is SRH's linearization, or 0 without SRH) and L_eff = 1/κ:
  - **Y_fixed = q D n₀ e^v/V_T · κ coth(κW′):** Shockley eq 4.23 with the finite base (PA; the coth form PD). It is G₀√(1 + iωτ) for a long base and the short-base series for a short one.
  - **Y_edge = q D U₀ (dW′/dV) [u₀″/U₀ − κ coth(κW′) u₀′/U₀]:** the depletion edges' motion, with dW′/dV = C_j/(qN). This is derived here (PD) and is not in eq 4.23. At ω = 0, Y_fixed + Y_edge is d J_diff/dV of diode.hpp (its W(V) edges) exactly.
  - **iω C_j,** with qe_junction's C_j at V.
  - The **Shockley total** (ΣY_fixed + iωC_j) and the **total** (adding ΣY_edge) are both reported.
  - Complex double-double is written locally, with cos and sin reduced by 2π and halved, an expm1 without cancellation, and coth for Re z > 0.
  - **Second method** for Y_fixed: the Riccati form w = δu/δu′, w′ = 1 − κ²w², integrated by Dormand–Prince in double-double from the contact (w = 0) to the edge, in y/|L_eff| over at most 20 |L_eff| (the contact's influence is e^−40). The difference is reported.
  - **A departure from 16.3, made in the build:** 16.3 names complex Chebyshev collocation as the second method. At high ωτ the base holds up to about 100 |L_eff|. Collocation would then need about 10 patches and a dense double-double solve of about 500 unknowns per frequency (about 1 s each), or a block-banded solver that does not exist. The Riccati form is a different formulation of the same ODE, independent of the closed form, and costs milliseconds. Your call if collocation is wanted anyway.
- **L0 (`admittance_test`, 7 cases, 151 assertions):**
  - at ω = 0, the total against a Richardson-extrapolated double-double difference of diode_current's J_diff, to 1e-9: with and without SRH, and for 5 µm and 500 µm bases. The edge term is resolved, between 1e-9 and 1e-2 of G;
  - the long-base limit at W′ = 40 L, to 1e-13;
  - the short-base series 1/W′ + κ²W′/3 − κ⁴W′³/45;
  - the Riccati method against the closed form over ωτ ≈ 8e-4 … 8e3, to 1e-12;
  - the edge term at ω > 0, against a second algebraic form built from Y_fixed;
  - iωC_j and the sums; the refusals.
- **Mutation:** 24 mutants, 22 caught on the first run. The survivors:
  - M22 (the sign of the edge term's imaginary part, which was checked only at ω = 0) is caught after adding the ω > 0 case;
  - M10 (the square root's branch for Re z < 0) was unreachable, since Re κ² ≥ 0. It is removed, and the precondition is stated.
- **Timing:**
  - Release: the 7 cases take 0.36 s.
  - Debug: `verify.admittance` takes 1.45 s. I trimmed the SRH biases in the ω = 0 case and removed duplicate operating-point calls to get it there.
  - The serial Debug verify run measured 95.4 s (66/66 pass), but with a game using about 56 % of the CPU: every unchanged test was 20–35 % slower than in 16.27 (fermi_dirac 21.3 → 25.4 s). The idle estimate is 72.5 + 1.45 ≈ 74.0 s against CV-Q19's 75 s, to be remeasured on an idle machine before close-out (step 14).

**Measured sizes for TVB-18's construction** (the reference itself, NiTCAD's spec values as TVB-16, τ_p = 8.33 µs the p base's minority lifetime, R = ΣW′/(qNμ_maj) with Caughey–Thomas μ_maj 405 and 1157 cm²/Vs). This was a temporary diagnostic, not committed.

| base | V | SNS/J_diff | \|Y\|R at ωτ = 1e-3 / 1 / 1e3 / 1e4 | \|Y_edge\|/\|Y\| at low ω | iωC_j share of \|Y\| at ωτ = 1 |
|---|---|---|---|---|---|
| 5 µm | 0.35 | 6.3e-3 | 2.9e-6 / 5.2e-6 / 4.3e-3 / 4.3e-2 | 1.3e-3 | 0.82 |
| 5 µm | 0.40 | 2.7e-3 | 2.0e-5 / 2.1e-5 / 4.6e-3 / 4.6e-2 | 1.4e-3 | 0.22 |
| 5 µm | 0.45 | 1.2e-3 | 1.4e-4 / 1.4e-4 / 5.1e-3 / 5.1e-2 | 1.5e-3 | 0.036 |
| 500 µm | 0.50 | 1.3e-2 | 4.0e-3 / 4.6e-3 / 0.65 / 5.9 | 3e-7 | 0.12 |
| 500 µm | 0.55 | 5.9e-3 | 2.8e-2 / 3.1e-2 / 1.3 / 8.1 | 3e-7 | 0.021 |
| 500 µm | 0.60 | 2.6e-3 | 0.19 / 0.21 / 5.5 / 21 | 4e-7 | 3.3e-3 |

- The neutral bases' series resistance enters Y as about |Y|R, and at high ω it is mostly ωC_jR. Neither Y_ref nor 16.18 (g)'s ε_Y contains it (finding (x)'s cause, now in AC).
- In the 500 µm base, every bias that meets 16.2 (l)'s SNS ≤ 1e-2 has |Y|R ≥ 4e-3 at low ω, and above 0.6 at ωτ = 1e3. The product (SNS)² × G₀R does not depend on V, so no bias satisfies both.
- In the 5 µm base, |Y|R ≤ 1.4e-4 up to ωτ = 1 and ≈ 5e-3 at ωτ = 1e3, against ε_Y's 3L_D/|L_eff| = 2.4e-2 there.

**TVB-18's construction, planned (my choices; nothing has run)**, assuming DT-Q29 A:
- The 5 µm device of TVB-16 (ii) short (SRH on, uniform nodes, ladder 1024 … 128, the factorial at 2048), at forward biases 0.35, 0.40 and 0.45 V.
- ωτ_p ∈ {1e-3, 1e-2, 1e-1, 1, 10, 1e2, 1e3}, from NiTCAD's `solve_small_signal` (f = ω/2π). Y_sim is Y₁₁ of the anode, A/(V cm²).
- |Y| and arg Y are judged by 16.19's rule, with ε_Y from 16.18 (g) measured on the solved state (Δn/N, the DT-Q25 ε_φ, J_SNS/J_diff, 3L_D/|L_eff|) and applicability |L_eff| ≥ 30 L_D (≥ 4.2 here) and ε_Y ≤ 0.1|Y_ref|.
- U_ref is the Riccati difference plus C_j's complex-step difference. U_num follows DT-Q26 and DT-Q28.
- 42 comparisons (3 biases × 7 frequencies × magnitude and phase); with DT-Q31 A, 63.

**DT-Q29 (the 500 µm base cannot meet 16.18 (g)).**
- (A) Build TVB-18 on the 5 µm base only, at 0.35–0.45 V. Record the long base as not buildable under 16.18 (g) on NiTCAD's silicon, with the table above. The √(1 + iωτ) behaviour is then seen only partly: at ωτ = 1e3, |L_eff| ≈ W′, the transition from short to long base.
- (B) Keep the long base and add the bases' impedance to Y_ref: Y = 1/(1/Y_j + Z_R), with Z_R from the reference's own minority profile (the majority current's ohmic drop). This is new reference physics and changes 16.3 and 16.18 (g) after their freeze (before any run). It is first-order only where |Y|R ≪ 1: it neglects the minority carriers' drift in the ohmic field, which is of the same order. So the long base above ωτ ≈ 1 would still be outside it.
- (C) Keep the frozen plan and run it. The long base fails like TVB-16's (x), and the fails become findings.

**Recommendation: A.** The table shows the long base fails for a reason 16.18 (g) does not budget, before anything is run. A avoids building a known finding, and B's extension would not hold where the long base matters.

**DT-Q30 (the edge motion in Y_ref).**
- (A) Y_ref = the total, with Y_edge. Its ω = 0 limit is d J_ref/dV of TVB-16's own reference, which uses W(V) edges (DT-Q5), and it is checked at L0. Shockley's fixed-edge form is printed beside it.
- (B) Y_ref = eq 4.23 with fixed edges, 16.3's literal form. Y_edge (1.3–1.5e-3 of |Y| in the 5 µm base) is then an unnamed difference inside ε_Y.

**Recommendation: A.** It is the exact linearization of the reference's own DC model and adds no new approximation. I-9's "eq 4.23's full form" (both carriers, plus C_T) is kept, and extended consistently with the W(V) edges.

**DT-Q31 (what magnitude and phase test at high ω).** In the 5 µm base at 0.35–0.40 V, iωC_j is 0.2–0.8 of |Y| at ωτ = 1, and nearly all of it at 1e3. There, magnitude and phase mostly re-test C_j (TVB-17's quantity) and see the diffusion admittance only weakly.
- (A) Also judge Re Y, the diffusion conductance G(ω), which C_j does not enter, by the same rule: |Re Y_sim − Re Y_ref| ≤ U_num + U_ref + ε_Y. This adds 21 comparisons and widens nothing.
- (B) Magnitude and phase only, as 16.4 states.

**Recommendation: A.** It makes the benchmark test what it is named for, at all frequencies.

**16.35 Rulings DT-Q28 B (authorized), DT-Q29 A, DT-Q30 A, DT-Q31 A** (owner, "choose recommended" 2026-11-06, and "Authorize B" on DT-Q28, asked separately because B widens U_num after TVB-16's run, TVM2)**; DT-Q28 B applied; TVB-16's frozen outcome.**
- **DT-Q28: B, authorized.** In `tvb16.cpp`, `u_disc()` is the larger of:
  - the in-range three-grid GCI (DT-Q26 A);
  - Roache's two-grid estimate of the N = 1024 value, F_s |f(h/2) − f(h)| r^p/(r^p − 1) with F_s = 3 (`gci_safety_factor_two_grids`) and p = p_th = 2 on N = 2048 and 1024 (the factorial's (h/2, tol) state). It is used only where the measured e(h/2)/e(h), the errors against the reference, is within 0.25 ± 0.05.

  If neither applies, U_disc is not estimated, and the reason names both causes. This holds for J and for the derivatives. The note prints which estimate was used.
- **DT-Q29: A.** TVB-18 is built on the 5 µm base only. The 500 µm base is recorded as not buildable under 16.18 (g) on NiTCAD's silicon (16.34's table).
- **DT-Q30: A.** Y_ref is the total, with Y_edge. Shockley's fixed-edge form is printed beside it.
- **DT-Q31: A.** TVB-18 also judges Re Y, by the same rule. That makes 63 comparisons.
- **TVB-16 rerun** (Release, first run after DT-Q28 B, nothing changed after it): **13 pass, 28 not_evaluated, 0 fail.** `verify_benchmark.tvr4_full` passes, with the frozen outcome set to this result: the 13 passing comparisons are listed in `device_theory_test`, and every other one must be not_evaluated.
  - **Pass:**
    - (i) J at 0.2 V (three-grid in range);
    - (ii) short: J at 0.2, 0.3, −0.5, −1 and −2 V, and the −0.5..−1 and −1..−2 V reverse exponents;
    - (ii) long: J at −0.5, −1 and −2 V, and the same two exponents.
    - The −2 V long-base J, 16.33's fail, now passes with the two-grid estimate (5.1e-4 against an error of 1.7e-4; 16.33's three-grid 1.15e-4 is the smaller of the two).
  - **not_evaluated, by cause:**
    - Most have e(h/2)/e(h) ≈ 1 (for example 1.0002 at 0.6 V and 1.011 at 0.5 V in (ii) short). There, the error against the reference is the reference's approximation (ε_J's terms, and in the long base finding (x)'s ohmic drop), not mesh error. The order test can't be applied and the fallback correctly declines.
    - Others have ratios far from 1/4 (−0.79, 2.19, 0.50, 1.66: errors at rounding or resolution level that change sign).
    - At −5 V, the factorial's tol/100 states were not reached (DT-Q23 A).
    - Some derivatives are not_evaluated because one of their J values is.
  - Finding (x) still classifies nothing: no long-base forward point is judged.
- **What TVB-16 now verifies:** reverse-bias J and its voltage exponent in both bases with SRH, and the low-forward J in the short base. In forward bias in the short bases, NiTCAD's J agrees with the reference within the bounds (log errors 2.4e-6 to 6.4e-3 against 1.3e-2 to 2.1e-2), but the discretization error can't be separated from the reference's approximation, so those points stay not_evaluated. The long base's forward J at 0.5 and 0.6 V would exceed its bounds (1.3e-3 and 5.7e-2 against 3.4e-4 and 7.3e-3: finding (x)'s ohmic drop), and is not_evaluated for the same reason. This is recorded as TVB-16's reach, not as a finding.



**16.36 Step 7b as built: TVB-18; the shared U_disc rule; TVB-16's mutation run. DT-Q32 and DT-Q33 open.**
- **The shared U_disc rule** (`registry.hpp`'s `u_disc`, beside `factorial`): the in-range three-grid GCI (DT-Q26 A), Roache's two-grid estimate gated as DT-Q28 B, and the larger of the two. TVB-16's local copy is replaced by it, and TVB-18 uses it. TVB-16's outcomes and its 13 pinned bounds are identical after the move. There are 8 unit cases in `registry_test` on exact sequences c + a h^p: values, orders, each reason, the band's edges, and the spread between 0.1 and 1.
- **TVB-16's mutation run** ("TVB-16*", 14 mutants): 8 caught. The 6 survivors each change U_num or ε but no outcome, because the 13 passes have wide margins:
  - M1, assumed order 4;
  - M3, F_s = 1;
  - M7, the coarser-triple check removed;
  - M8, h/2 halved;
  - M11, ε_φ dropped;
  - M12, the GCI denominator's sign.

  Changes made:
  - the frozen test now pins the 13 passes' bounds to 1e-6;
  - the rule moved to the registry, with unit tests.

  Rerun, the 6 survivors and 7 mutants of the shared rule: all caught. The spread and band survivors were caught after the two added cases.
- **TVB-18** (`tvb18.cpp`, registered after TVB-17, full tier, 47 benchmarks), as planned in 16.34 under DT-Q29–Q31 A:
  - 5 µm bases, the ladder 1024 … 128 and the factorial at 2048. Biases 0.35, 0.40 and 0.45 V. ωτ_p = 1e-3 … 1e3 with τ_p = 8.33 µs.
  - |Y|, Re Y and arg Y of Y₁₁, 63 comparisons.
  - The phase bound is asin((U_num,|Y| + U_ref + ε_Y)/|Y_ref|): the largest angle a complex error of the bound's size can produce (16.19's bound read on the complex value).
  - The frozen expected outcome, set before the first run, is "every comparison passes". The case is in `device_theory_test` [.tvr4_full], and the run takes 0.5 s in Release.
- **TVB-18's first run** (nothing changed after it): **1 pass, 62 not_evaluated.** `tvr4_full` fails on TVB-18.
  - U_disc is not estimated at 62 points:
    - three-grid outside the range at 37, no coarser order at 18, oscillatory at 7;
    - the two-grid gate is refused everywhere, with e(h/2)/e(h) = 0.97–1.0.
  - The level-to-level differences are 1e-6 to 4e-4 of |Y|. The error against the reference (1e-3 to 3e-2) is the reference's approximation and NiTCAD's model departures, not mesh error. So the gate's ratio is about 1 wherever ε_Y dominates.
  - This is a flaw in the gate I proposed in DT-Q28 B, not in NiTCAD. It is also why TVB-16's forward points are not_evaluated (16.35).
- **Two physical observations from the same run:**
  - **(a) Series resistance:** at ωτ = 1e3, Re Y_sim / Re Y_ref = 7.0, 2.0 and 1.19 at 0.35, 0.40 and 0.45 V. The excess matches ω²C_j²R (R = 8.8e-4 Ω·cm², 16.34). This is the term DT-Q29 kept out of reach. Re Y's bound is ε_Y|Y_ref| in absolute terms, so it covers it.
  - **(b) C_j above the QE reference:** Im Y_sim / Im Y_ref = 1.0068, 1.0133 and 1.0272 at 0.35, 0.40 and 0.45 V, the same at ωτ = 1e-3 and 1e3 (Im Y ≈ ωC_j throughout, since the short base's diffusion susceptance is about 1e-3 of it).
    - The excess roughly doubles every 50 mV and is converged in the mesh (e(h/2)/e(h) ≈ 1). At reverse bias TVB-17's C agrees to about 1e-5.
    - None of ε_Y's named terms has this dependence. Δφ/V_T is 4.7e-3, flat in V. J_SNS/J_diff falls with V.
    - Not explained yet.
- **Diagnostic** (a temporary build, reverted): the observed order on the finest triple N = 2048/1024/512 (the factorial's (h/2, tol) state with levels 0 and 1).
  - In TVB-18 the median order is 2.04 (range −0.9 to 2.41), and 7 of 63 change sign.
  - With the two-grid estimate gated on that order instead (within 2 ± 0.5, differences of one sign):
    - **TVB-18 would be 54 pass, 1 fail, 8 not_evaluated.** The fail is |Y| at 0.45 V, ωτ = 100: error 1.53e-2 against a bound of 1.29e-2, from (b). The 8 are |Y| and Re Y at 0.35 V with ωτ ≤ 10, where the differences change sign.
    - **TVB-16's J would be 14 pass, 2 fail, 7 not_evaluated.** The 2 fails are the long base at 0.5 and 0.6 V, inside finding (x)'s region, so they become findings as DT-Q27 intended. TVB-16's derivatives were not measured.

**DT-Q32 (the two-grid gate).**
- (A) **Gate the two-grid estimate on the discretization itself:** the observed order on N = 2048/1024/512 within p_th ± 0.5, with successive differences of one sign, instead of e(h/2)/e(h) against the reference. The rest of DT-Q28 B (F_s = 3, p_th = 2, the larger of the two, not estimated where neither) is unchanged.
  - Applied to TVB-16 and TVB-18. Both have run, so this is a post-run change (TVM2) and needs your **authorize**.
  - The ±0.5 band is mine. 4.2's 0.1 spread is the three-grid rule, which these ladders don't meet.
- (B) **Keep DT-Q28 B as ruled.** TVB-18 is frozen as found (1 pass / 62 not_evaluated), and TVB-16 stays at 16.35.
- (C) **The two-grid estimate with F_s = 3 and no gate.** That gives TVB-18 62 pass / 1 fail, but no check that p_th = 2 holds (the 7 sign changes included).

**Recommendation: A.** It is what DT-Q28 B meant to test, whether the mesh error falls as h², measured where it can be. B leaves both benchmarks verifying almost nothing, and C trusts p_th without looking.

**DT-Q33 (C_j above the QE reference at forward bias, 0.7–2.7 %).**
- (A) **Investigate before TVB-18's frozen outcome is set.** This is a diagnostic only, with no criterion change:
  - NiTCAD's quasi-static capacitance −dQ/dV from its own DC states (the charge in the depletion region and in each neutral base, separately), against qe_junction's Q(V) and the reference's diffusion charge;
  - whether the excess sits in the depletion region (the QE model's charge) or in the bases (charge the reference does not have).

  I report the result and propose either a reference correction (as its own ruling) or a finding.
- (B) Record it now as a known finding on TVB-18's |Y| where it decides the outcome (one point under DT-Q32 A), cause unknown.

**Recommendation: A.** A 2.7 % capacitance difference that grows with bias is either a reference term missing at forward bias or a NiTCAD defect. Both matter commercially, and the two are distinguishable from charges already on hand.


**16.37 DT-Q32 A and DT-Q33 A (authorized: owner "go ahead", 2026-11-06, in reply to "authorize DT-Q32 A, DT-Q33 A"); TVB-16 re-frozen; the C_j excess explained. DT-Q34 open.**
- **DT-Q32 A applied** (a post-run change, authorized as above). The shared `u_disc` now gates the two-grid estimate on the discretization itself: the observed order on N = 2048/1024/512 (the factorial's (h/2, tol) state and ladder levels 0, 1), within 2 +- 0.5, with differences of one sign. The rest of DT-Q28 B is unchanged (F_s = 3, p_th = 2, the larger of the two, nothing where neither is kept). `u_disc` takes no reference errors any more; `Discretization` gains `order_fine`.
  - `registry_test`'s U_disc case is rewritten for the gate: exact order 2 gives the two-grid estimate; finest-triple orders 1 and 3 are not kept; an oscillating finest triple; band edges 1.45/1.55/2.45/2.55; the spread cases. 50 assertions pass.
- **TVB-16 rerun and re-frozen as found:** 19 pass, 3 findings (x), 19 not_evaluated.
  - The 3 findings are J (ii) long base at 0.5 V and 0.6 V and the ideality between them: the series-resistance finding, as DT-Q27 A intended.
  - The not_evaluated points: the finest-triple order outside 2 +- 0.5 (0.65–4.9) or oscillating, with the three-grid GCI also out of range; or a factorial state not reached.
  - The frozen case pins the 19 bounds to 1e-6.
- **TVB-18 rerun:** 51 pass, 1 fail, 11 not_evaluated.
  - The fail is |Y| at 0.45 V, ωτ = 100: 1.53e-2 against a bound of 1.29e-2.
  - The not_evaluated points are all at 0.35 V: 8 |Y| and Re Y points whose differences change sign, and 3 arg Y points whose bound takes |Y|'s U_num.
  - The frozen case ("every comparison passes") is unchanged and still fails. It is not re-frozen pending DT-Q34. `tvr4_full` fails on TVB-18 only. The default verify suite passes, 68 of 68 outside `tvr4_full`.
- **DT-Q33 A, the investigation** (a temporary diagnostic, reverted; no criterion changed). Using NiTCAD's own quasi-static (f = 0) small-signal fields at N = 1024 and 2048, 0.1–0.45 V:
  - **The depletion charge agrees.** −dQ_n/dV, the n half's net charge response (depletion plus tail), equals qe_junction's C_j to about 1e-4. The excess is not in the QE charge.
  - **The base storage agrees.** The minority storage in each base, times the short base's 2/3, matches the reference's Im(Y_fixed + Y_edge)/ω.
  - **The cause.** The terminal current at x = 0 is exactly iωC_j + J_n(0) + J_p(0). The reference takes the minority currents at the depletion edges, x_p and x_n, not at x = 0. So it omits the minority carriers stored inside the depletion region, C_SCR = q(∫_0^x_n ∂p/∂V dx + ∫_x_p^0 ∂n/∂V dx). Its in-phase counterpart, the SNS recombination, is in ε_Y; this quadrature term is not.
  - **The measure.** At N = 2048, (Im Y_sim/ω − C_ref,total) / C_SCR is 1.021, 1.008, 1.002 and 0.9996 at 0.30, 0.35, 0.40 and 0.45 V. At lower bias the residual is a few 1e-12 F/cm², and it falls with h: 1.31 → 1.09 at 0.1 → 0.2 V for N = 2048, against 2.14 and 1.33 at N = 1024.
  - So the excess is a term the reference omits, not a NiTCAD defect. C_SCR / C_j is 0.64 %, 1.29 % and 2.70 % at 0.35, 0.40 and 0.45 V.

**DT-Q34 (the depletion region's minority storage, missing from TVB-18's reference).**
- (A) **Reference correction:** add iωC_SCR to the reference.
  - C_SCR is computed in the QE model, in double-double: p and n are referred to the flat quasi-Fermi levels, x(u) comes from the first integral by quadrature, and the integrals run between the charge-balance edges and the metallurgical junction.
  - The second method is complex step on the same integrals, with admittance_test cases.
  - A post-run change to the reference; needs your **authorize**. TVB-18 would then be frozen as found after one rerun.
- (B) Add C_SCR/C_j to ε_Y as a named approximation term (a post-run widening; needs **authorize**).
- (C) Keep the reference; record a known finding on TVB-18 |Y| at forward bias where C_SCR decides the outcome (the one point), with this cause.

**Recommendation: A.** The term is physical and computable in the same QE model. Adding it tightens the check of NiTCAD instead of widening a bound, and the diagnostic shows NiTCAD agrees with it to under 1 % at 0.35–0.45 V.

**Not done yet** (after DT-Q34): freeze TVB-18 and run its mutants ("TVB-18*"); the mutants of the new gate ([registry]) and "TVB-16*"; then step 8 (TVB-10).

**16.38 DT-Q34 A (authorized: owner "go A", 2026-11-06): the reference gains i w C_SCR; TVB-18 frozen as found.**
- **The correction** (a post-run change to the reference, authorized). `verify/reference/admittance` adds i w C_SCR to Y_ref's total. C_SCR = q(∫_0^{x_n} ∂p/∂V|_x dx + ∫_{−x_p}^0 ∂n/∂V|_x dx) over the QE profile at V, between the charge-balance edges and the junction. Shockley's form is reported without it. The new fields are `DiodeAdmittance::scr_storage_F_per_cm2` and `scr_second_method_difference`, and `AdmittancePoint::scr_storage`.
  - **Derivation** (PD; `admittance.hpp`). At fixed x, ∂u/∂v = −(∂x/∂v|_u)·√(kΦ). Exchanging the order of the double integral leaves single integrals in u, and the inner minority integral is closed-form. The closed forms are ∂Φ_b/∂v|_u = (n + p − n_b − p_b)/2 (from Φ_b = n + p − (n_b + p_b) − C(u − u_b)) and du_j/dv = (√(N_A² + 4m²) − √(N_D² + 4m²)) / (2(N_A + N_D)). The integrals use Gauss-Kronrod in u, in double-double.
  - **Departure, second method.** DT-Q34 A named complex step on the same integrals. The QE profile (`pn_potential`'s Newton inversion and Brent's root) has no complex extension here, so the second method is a different one: the minority charge between the fixed positions −x_p and x_n (at V), evaluated as a function of the bias and differenced by Richardson's 4-point rule in v (h = 1e-3). It shares the profile and quadrature, but not the differentiation or the change of variables. The difference is reported.
- **L0** (`admittance_test`, a new case). At 0.35 and 0.45 V, C_SCR is positive, grows with bias, and stays below 5 % of C_j. Its second method agrees to within 1e-9 relative at 0.45 V. At −2 V, C_SCR is below 1e-9 C_j and the second method agrees. The totals case now includes the term. 8 cases, 167 assertions pass.
  - **Cost (CV-Q19).** In the first build every `diode_admittance` call also ran the second method. That took `verify.admittance` from 1.45 to 5.2 s, and serial Debug verify to 80.6 s, against the 75 s budget.
    - The second method is now its own function, `scr_storage_second_method`, called by the L0 test and by TVB-18's notes. `diode_admittance` keeps only the first method, with its two integrals per side merged into one quadrature.
    - The new case is trimmed to 2 forward biases, and the quadratures run at 1e-20. That is still far below any criterion: C_SCR is a few percent of C_j.
    - Now `verify.admittance` takes 2.1 s, and serial Debug verify passes 66/66 in 74.2 and 74.6 s.
  - Mutation, on the final code: 8 mutants (the ½ factor, the u_j-term sign, the n and p sides' b-term each dropped, du_j/dv swapped, ∂Φ/∂v's ½, the second method's difference rule, the total without the term). All 8 are caught.
- **TVB-18's notes** report C_SCR per bias: 0.632 %, 1.28 % and 2.69 % of C_j at 0.35, 0.40 and 0.45 V. This agrees with 16.37's diagnostic (0.64, 1.29 and 2.70 %). The second method differs by about 1e-23 F/cm², 1e-14 relative.
- **Suites.** Release verify passes 74/74, the local full tier included (2622 s; `verify_full.tvb04` 1272 s); this ran before the restructure. After the restructure: `[admittance]`, `[tvr4_full]` and serial Debug verify pass again.
- **TVB-18 rerun: 52 pass, 11 not_evaluated, 0 fail** (it was 51/11/1). The |Y| miss at 0.45 V, ωτ = 100, is gone. The 11 not_evaluated points are unchanged: at 0.35 V, |Y|, Re Y and arg Y at ωτ_p = 0.001–0.1, and Re Y at 1 and 10, where U_num is not estimated.
  - **Frozen as found** (DT-Q34 A): `device_theory_test`'s TVB-18 case pins those 11 as not_evaluated and the other 52 as pass. The criteria are unchanged; only the reference moved.
  - `tvr4_full` passes: 6 cases, 170 assertions (Release).
- **Not done yet:** the mutants of TVB-18's own code ("TVB-18*"), of the DT-Q32 gate ([registry]) and of "TVB-16*"; then step 8 (TVB-10).

**16.39 Step 7's mutation runs (owner "complete TVR-4 in one go", 2026-11-06).** No criterion, reference or outcome
changed. Each mutant was applied to the Release build, the frozen case run, and the file restored.
- **TVB-18's frozen case now pins its 52 passes' bounds to 1e-6**, as TVB-16's does (16.36). The outcomes were already
  pinned; the bounds were not, so U_num and ε_Y mutants could only be caught through an outcome change. This is test
  strengthening only (the values are those of 16.38's run).
- **TVB-18\*** (17 mutants: both applicability constants, each ε_Y term dropped or rescaled, L_eff's min, U_ref's
  complex-step term, the phase bound's U_num, the factorial's states, the resolution, the edge densities' doping,
  the quasi-Fermi factor, both h's, Y₁₂ for Y₁₁, the judge's bound): **13 caught.** The 4 survivors are equivalent on
  this data:
  - U_ref's ω|C_j − C_j,cs| term (U_ref is 3e-18 S/cm², below the pins' 1e-6);
  - the factorial's (h/2, tol/100) state for (h/2, tol), and (h, tol) for (h, tol/100): Newton converges past both
    tolerances, and the states are bit-identical (ε_int = 0 at every point);
  - the bound halved: the largest metric/bound is 0.143.
- **[registry], the DT-Q32 gate** (11 mutants: the band, its `<=`, the two-grid order and safety factor, the triple's
  h and f, the flat-triple clause, the oscillation test, max for min, the u2-only branch, the two-grid difference):
  **9 caught, then 10.** The flat-triple clause (`*c == 0.0`) survived; a `registry_test` case was added (finest
  triple exactly flat, the ladder's GCI refused: the two-grid estimate is kept and is 0), and the rerun caught it. The
  `<=` → `<` mutant is equivalent: the iterated order never lands exactly on 1.5 or 2.5.
- **TVB-16\*** under the DT-Q32 gate (14 mutants: h/2, each ε_J term, the edge factor, the resolution, U_ref's
  quadrature term, U_num in the bound, the derivatives' ε and resolution, the factorial's h/2 state, h, the
  derivative share, the applicability predicate, the factorial's argument order): **11 caught.** One mutant (the
  quasi-Fermi factor to 1) did not compile (unused variable, warnings as errors); its replacement was caught. The 3
  survivors are equivalent:
  - U_ref's quadrature term (U_ref/J is the 1e-15 floor at every bias);
  - the predicate forced true: ε_J ≤ 0.1 J_ref holds at every point of this data;
  - (h, tol/100) swapped with (h/2, tol): ε_int is symmetric in them.
- Release `[registry]` (53 assertions) and the TVB-16/18 frozen cases pass after the restores.

**16.40 Step 8 as built: mos_hf and TVB-10. DT-Q35 open.**
- **`verify/reference/mos_hf`** (PD, base-free; p-type substrates only, an n-type one is refused):
  - The redistribution model of 16.18 (b): the electrons at one quasi-Fermi level shifted by v, Φ(t, v) =
    p_b G(−t) + n_b e^{−v} G(t). This is 16.17's reproduction (2) with χ = 0. Its double zero at t = 0 keeps
    Q_n = −q ∫ n_b e^{−v} expm1(s)/g ds finite. Reproduction (1) is not used (16.17's finding).
  - C_s,HF = −(1/V_T)(∂Q_s/∂t + ∂Q_s/∂v · dv/dt), with dv/dt = −(∂Q_n/∂t)/(∂Q_n/∂v). The partials are analytic, and
    ∂Q_n/∂v is a Gauss–Kronrod quadrature in dd.
  - **Departure from 16.18 (b):** the plan named complex step for the partials. They are analytic instead, so the
    second method (Q_s along Q_n held fixed by Brent, Richardson central differences, δ = 1e-3) is independent of
    them.
  - 16.18 (b)'s "V frozen" diagnostic is C_LF itself: v held at 0 means the electrons follow the substrate. The
    benchmark prints C_LF beside C_HF instead.
  - Deep depletion: Φ_dd = N_A G(−t), with V_G by Brent. The second method is the quadrature of the charge density
    over the profile.
- **L0** (`mos_hf_test`, 5 cases, 71 assertions, 0.25 s Debug):
  - v = 0 is KN's charge (1e-28);
  - Q_s = Q_p + Q_n by quadrature, inversion and accumulation, v ≠ 0 (1e-24);
  - depletion C_HF = C_LF (1e-6); the strong-inversion plateau within 1 %, and inside ±10 % of ε/W_max;
  - deep depletion: its second method (1e-24), the round trip (1e-25), KN in depletion (1e-8), and smaller |Q|
    with larger band bending in strong inversion;
  - the refusals.
  - The second method takes 0.7 s in Debug against CV-Q19's exhausted budget. It is therefore not run in L0.
    TVB-10 runs it at all six gate biases and stops with a run error unless it agrees to 1e-15 C. The agreement
    found is 1e-17 to 1e-28 relative (6.8e-28 F/cm² at 2 V).
- **TVB-10** (`tvb10.cpp`, full tier, registered after TVB-09; 48 benchmarks). TVB-09's MOS-C, now with SRH on,
  Caughey–Thomas mobility and no field mobility.
  - HF (A): the operating point is `equilibrium_poisson`, as `small_signal.hpp` prescribes for inversion. The small
    signal runs at 1e5, 10^5.5 and 1e6 Hz, and C = Im Y_gg/ω is judged at 10^5.5 Hz. The ladder is N = 512…64, with
    N = 1024 as h/2. U_disc comes from the shared `u_disc`, and the (h, tol) factorial runs at N = 512.
  - Deep depletion (E): a gate ramp from V_FB to V_FB + 4 V in 100 µs, BDF2 with fixed 1 µs steps, at Newton
    tolerance 1e-10 (1e-12 for the factorial's tol/100 state). Q_s is judged at 1, 2, 3 and 4 V. The factorials are
    (h, tol) and (h, Δt) at N = 512, with Δt/2 the finer step.
- **Construction changes before commit** (no criterion, ladder or expected outcome touched; as TVB-15's, 16.26):
  - The first build used a 1 µs ramp with adaptive steps (rtol 1e-6, factorial rtol 1e-8). The factorial's
    adaptive runs took more than 45 minutes and were killed.
  - Fixed steps of 10 ns made the run 2 s. That run found the holes' lag behind the ramp: Q_s at a ramp 10× slower
    differed by 1.7e-6. That is comparable to the finest mesh error (3.7e-6), and it bent the 1 V order to 2.46.
  - The ramp became 100 µs. That is still deep depletion: generation needs about 0.1 s, and the electrons' charge
    measured at 4 V is 6.5e-7 of |Q_s|.
  - The Newton tolerance became 1e-10, because tol/100 = 1e-14 fell below the rounding floor at 1 µs steps
    (Newton failed at t = 41 µs).
  - The deep-depletion truncation uses the bound 4 e^{−(d − W)/L_D} (at most 3e-10). `bulk_profile` cannot take a
    deep-depletion band bending: it would include the inversion layer.
- **Frozen outcome, set before the first run:** every comparison passes (14 comparisons).
- **TVB-10's run: 13 pass, 1 not_evaluated, 0 fail.** 4.3 s in Release. The frozen case fails, so `tvr4_full`
  fails, pending DT-Q35.
  - HF C_HF passes at V_G − V_FB = 0.2, 0.5, 1, 1.5 and 2.5 V, with |C_sim − C_ref| = 4.3e-13 to 1.0e-12 F/cm²
    (5.6e-6 to 1.2e-5 relative). ε_HF, the spread over the decade, is 9e-18 to 1.2e-12, far inside 1e-3 C.
  - **2 V is not_evaluated.** The three-grid GCI is outside the asymptotic range (orders 2.108, 2.005), and the
    finest triple oscillates. The (h/2, tol) state is not converged in C: e(h/2, tol) = 1.43e-11 against
    e(h/2, tol/100) = 4.26e-12.
  - Deep depletion: all 8 comparisons pass. p̂ = 2.001, 1.975, 1.941 and 1.903 (spread ≤ 0.071). The finest errors
    are 5.6e-6 to 6.4e-6. U_iter is ≤ 1.2e-12, U_dt ≤ 5.1e-11 and ε_int ≤ 1.6e-12. The orders' drift toward 1.9 at
    high V_G follows the electrons' growing charge (1.6e-8 to 6.5e-7); the slower ramp's departure is 7.8e-8 to
    4.8e-6.
  - **Observation (not a finding, no criterion involved):** in strong inversion NiTCAD's reported small-signal
    resolution is 3.2e-10 to 4.7e-9 F/cm² (up to 6 % of C). The actual agreement with the reference is 5e-13. The
    resolution dominates U_num there, so the HF bounds at 1.5–2.5 V are 100–1000× wider than the error.
- **Mutation:**
  - `mos_hf` (10 mutants: Φ's shift, Φ_v's form, the ½ in Φ_v/2Φ, dv/dt's sign, Richardson, the deep-depletion
    charge and density, dQ_n/dt's sign, the v-term): 9 caught.
    - The Φ_v mutant first survived the 1e-9 second-method bound. The bound was tightened to 1e-15 (the agreement
      found is ~1e-17), and it was then caught.
    - The survivor (dQ_n/dt's sign, which acts only in accumulation) is equivalent: the term it re-signs is about
      1e-17 of C there.
    - After the second method moved into TVB-10, the Φ_v and Richardson mutants were rerun against TVB-10. Its
      self-check catches both.
  - TVB-10 (14 mutants): judged by any change in its 14 printed comparison lines, because the frozen case is
    already red. 12 caught. The 2 survivors are equivalent:
    - U_dt dropped: it enters only the factorial's interaction test, which is far from its threshold;
    - the factorial's argument order: ε_int is symmetric in it.

**DT-Q35 (TVB-10's frozen outcome).**
- (A) Record as found, as DT-Q20 A did: the 2 V HF point is not_evaluated with its reason (the h/2 state's tolerance
  sensitivity), and the other 13 pass. The frozen case then pins those outcomes and the passes' bounds to 1e-6.
  Criteria, ladder and construction are unchanged.
- (B) Tighten the HF Newton tolerance to 1e-14 (factorial 1e-16, if reachable) and rerun. That changes the
  construction after a run (TVM2), and 1e-16 is probably below the rounding floor.

**Recommendation: A.** The miss is the U_num estimate at one point, not an error: that point's |C_sim − C_ref| is
5.3e-13, the same as its neighbours.

**16.41 Step 10 as built: linear_transient and TVB-19. DT-Q36 open.** (Step 9, Pao–Sah, is taken after the lighter
steps; 16.10's order is otherwise kept.)
- **`verify/reference/linear_transient`** (PD, base-free):
  - `dielectric_relaxation_time`: τ_d = ε/(q μ N) in dd. With σ = q(μ_n n_0 + μ_p p_0), the benchmark passes μ = 1
    and N = σ/q.
  - `relaxation_second_method`: Dormand–Prince with 4000 fixed steps on dQ/dt = −Q/τ_d, against the exponential, to
    5 τ_d. The largest difference is 3.1e-19.
  - `haynes_shockley_moments`: the moments of a drifting, spreading, decaying Gaussian (σ, D, v, τ), integrated over
    the source's time profile. Its age a lies on [t − t_off, t − t_on] with the weight e^{−a/τ}, so the moments are
    those of a truncated exponential, in closed form:
    - area = G σ √(2π) e^{−a₁/τ} m₀;
    - centre = x₀ + v⟨a⟩;
    - width² = σ² + 2D⟨a⟩ + v² Var(a).
  - **Departure from 16.3:** the plan's second method was the PDE by method of lines with Dormand–Prince. It is
    instead Gauss–Kronrod on the three age integrals. A method-of-lines solution would carry its own spatial error,
    at the level being judged. The quadrature checks the closed forms independently: it agrees to 4.5e-33, 4.3e-31
    and 7.5e-29 at the benchmark's point. The Gaussian step (variances add under drift and diffusion) is checked by
    L0's limit cases, not by the second method.
- **L0** (`linear_transient_test`, 4 cases, 40 assertions, 0.16 s Debug): τ_d by hand; the second method; the
  uniform-age limit (τ = ∞); the moments' equations at τ finite against a direct quadrature; the short-pulse limit; the
  refusals.
- **TVB-19** (`tvb19.cpp`, registered after TVB-18; 49 benchmarks). The construction is as in the file header, from
  16.2 (m) and 16.8. Injection uses the TVQ3 hook's time profile (DT-Q3 A).
  - **Tier:** the full tier takes 438 s in Release: 16 runs per case, at N up to 4097, with fixed steps. That is past
    `tvr4_full`'s 120 s. It therefore runs in the local nightly tier, as `verify_full.tvb19` (DT-Q9: heavier ladders
    local). The fast tier has no comparisons.
- **Frozen outcome, set before the first run:** every comparison passes (16).
- **TVB-19's run: 10 pass, 3 not_evaluated, 3 fail.** The frozen nightly case fails, pending DT-Q36.
  - The finest accuracies all pass: decay rate 2.87e-5, area 5.6e-7, centre 4.1e-5 of the drift length, width
    8.5e-5, each against 1e-3. The backward-Euler orders pass (0.98–1.00). BDF2's dt orders pass for area and centre
    (2.000).
  - **Fail: the h order of the decay rate, area and centre** (p̂ = −0.0028, 3.2e-7, 2.0e-8). These quantities have
    no h error at the ladder's resolution. Their error is flat in h to 2e-8 of itself, and the error that exists is
    the time step's. The finite-volume sums conserve the integral, and the first moment under a uniform field. p̂
    is then undefined, and `judge_order` reports it as missed.
  - **Not evaluated: the decay rate's BDF2 dt order** (0.907 on the finest pair; 2.00 and 2.02 coarser). This is a
    floor of 2.86e-5 that neither h (2e-8) nor the tolerance (identical at tol/100) moves. Nor does dt below dt₀:
    the factorial's dt₀/2 state gives 2.92e-5.
    - A temporary diagnostic (not committed) varied the amplitude and N. At amplitude 1 and 0.25, N = 513 and 2049,
      rate·τ_d − 1 = 2.76e-5 and 2.87e-5, then 7.06e-5 and 7.46e-5. That fits c + d/amplitude, with c ≈ 1.4e-5 and
      d ≈ 1.4e-5.
    - The floor is not nonlinearity, which would grow with the amplitude, and not Newton.
    - The cause is not identified. It is reported as found: not a finding, and inside the accuracy criterion by
      35×.
  - **Not evaluated: the pulse width's BDF2 dt order (1.06) and h order (1.54).** Each ladder reaches the other
    parameter's floor. The width's finest error is mostly spatial (e(h/2) = 4.46e-5 against e(h) = 8.54e-5), so the
    dt ladder at the finest h flattens. The h ladder at dt₀ flattens on the time error.
  - The factorials are evaluated and non-interacting for all four quantities (|ε_int| ≤ 1.2e-7).
  - A cosmetic defect: the (h, dt) factorial's note keeps the shared formatter's labels "tol" and "tol/100" for dt₀
    and dt₀/2. It is left as is, to keep the run's output frozen.
- **Mutation:**
  - `linear_transient` (12 mutants: τ_d, the second method's rate, m₀, m₁ and Var's forms, the decay weight's age,
    the √(2π), the 2D, the centre's mean age, the uniform-age limit, the mean age's end, the refusal): 12 caught.
  - TVB-19: each mutant run takes about 7.5 minutes. It is judged by any change in the 16 printed comparison
    lines, because the frozen case is red. 3 mutants: the decay rate's τ_d scale (1e-4), the drift velocity's
    mobility (μ_p at 1e14 instead of 1e15) and the linearity predicate's threshold. All 3 caught.
  - The reduced set is a departure from the other steps' full mutation runs, made for cost. The judging code is
    shared with TVB-10's pattern, whose 14 mutants ran in 16.40.

**DT-Q36 (TVB-19's frozen outcome).**
- (A) Record as found, as DT-Q20 A and DT-Q35 A do: the frozen case pins the 10 passes (their bounds to 1e-6), the 3
  not_evaluated with their reasons, and the 3 h-order fails as fails with the reason above. Criteria, ladders and
  construction are unchanged.
- (B) Change the judging after the run (TVM2):
  - judge each one-parameter order on errors with the other parameter's part removed by Richardson;
  - judge the h-exact quantities as exactness (|e(h) − e(2h)| ≤ U_iter-based bound) instead of an order.
  This would turn the 3 fails into passes and probably 2 of the not_evaluated as well.

**Recommendation: A** for TVR-4. The fails are a property of the construction: the quantities chosen have no h error.
They are not errors in NiTCAD. B is the better design, but it changes criteria after seeing their results.

**16.42 Step 11 as built: energy_layer and TVB-13. DT-Q37 open.**
- **The linearization of 16.18 (a), checked before code** against ED2 (Stratton's current), Poisson and ED3 (the
  energy flux with c = 5/2, the Joule term J·E, the relaxation term): the three rows and the quartic P = A·D − B·C
  are as written there.
- **`verify/reference/energy_layer`** (PD, base-free):
  - `durand_kerner`: Weierstrass iteration on a monic real polynomial in complex dd, after scaling by the roots'
    geometric mean.
  - `energy_layer`: P's coefficients in dd, its four roots, λ = 2/min |Re s|.
  - **Departure from 16.18 (a):** the second method is not Chebyshev collocation of the linearized ODE with a decay
    fit. It is the determinant of the unreduced 3 × 3 system at each root, relative to the sum of its six expansion
    terms' magnitudes.
    - It is independent of the elimination and of the root finder. Collocation would have checked the same thing
      at much greater cost, and neither checks the linearization's signs against NiTCAD's equations.
    - That check is the printed diagnostic below: the decay measured on NiTCAD's own solution.
  - Found: 2e-32 to 1.7e-30 at the benchmark's three states. (The first normalization, by the rows' norms, gave
    1e-51, which is not a measure of cancellation. It was replaced before the first NiTCAD run.)
- **L0** (`energy_layer_test`, 4 cases, 37 assertions): a quartic with real and complex roots (1e-28); roots over ten
  decades (1e-26 relative); the zero-field limit ±1/L_D, ±√(G/κ) (1e-9); the determinant at three hot states
  (≤ 1e-25); the refusals.
- **TVB-13** (`tvb13.cpp`, Tier fast, registered after TVB-12; 50 benchmarks):
  - Drift (fast and full): field mobility on, L = 1 µm, N = 64, V/L = 1e2 … 1e5 V/cm.
  - Energy balance (full): N = 128 and 512, judged beyond 30 λ from each contact (13 ln 10 · λ).
  - μ₀ = the materials' μ_max (doping mobility off, as TVB-11). τ_w is silicon's cited 0.35 ps (ED5 as built: the
    default silicon set does not carry it). The reference values come from `model_mobility` in dd.
  - **Lengths chosen from λ before any NiTCAD run:** λ = 96 nm (1e3 V/cm), 0.20 µm (1e4) and 27 µm (1e5). The
    roots at 1e5 V/cm include 738 /cm, the unscreened cold-carrier mode that EC3's ruling (2026-10-29) describes.
    L = 20 µm for 1e3 and 1e4 V/cm, and 2 mm for 1e5.
  - **Construction change after a run error** (no outcome seen): at 1e5 V/cm the 2 mm bar stopped at 632 V (no
    convergence even at 2⁻²⁰ of the step). At the 20 kV the point needs, the potential's rounding (about 8e-11 in
    units of V_T) is above Newton's 1e-12 tolerance and comparable to the D criterion. 1e5 V/cm was removed from
    the energy-balance part; it remains in the drift part. 8 comparisons, not 10.
- **Frozen outcome, set before the first run:** every comparison passes.
- **TVB-13's run: 4 pass, 4 fail.** 0.06 s fast and about 1 s full, in Release. The frozen full case fails, so
  `tvr4_full` fails, pending DT-Q37.
  - Drift: all 4 pass, with relative errors 7.9e-14, 2.4e-14, 8.6e-16 and 1.7e-16 against 5.8e-11.
  - **Fail: T_n − T_L at both fields and both N.** The errors are 1.42e-4 and 7.3e-5 (1e3 V/cm, N = 128 and 512),
    then 6.6e-3 and 4.0e-3 (1e4 V/cm), against 2.3e-10 and 3.7e-9.
  - **Cause:** the reference was evaluated at V/L. With energy balance the contact layers carry part of the
    voltage, so the interior's uniform field is not V/L. Measured (printed, not judged): (E_int − V/L)/(V/L) is
    7.2e-5 and 3.7e-5 (1e3 V/cm), then 4.7e-3 and 2.9e-3 (1e4 V/cm). It falls with N because the layer is under-
    resolved at these meshes.
  - At each judged node's own field (its two edges' mean), the rise agrees with the homogeneous rise to 1.3e-14 …
    1.5e-13, inside 64 ε N². That is 16.2 (f)'s premise: the uniform state at its own field solves the discrete
    equations exactly.
  - This is a construction error in the benchmark (V/L was this step's choice; 16.4 says "T_n(E)"). It is not a
    finding.
  - **λ checked on NiTCAD's solution** (printed, not judged): |rise − rise(centre)| decays near the contact at
    2.03e5 and 9.5e4 /cm at N = 512, against min |Re s| = 2.09e5 and 9.9e4 /cm. At N = 128 it decays at 1.55e5
    and 7.2e4 /cm: under-resolved, but inside the ×2 safety and the 1e-13 margin. This confirms the linearization's
    signs.
- **Time:** the fast tier takes 1.97–2.59 s in Debug (0.06 s Release). With TVB-06's 0.11 s that is above 16.7's
  2 s for both cases. One solve per field instead of the continuation was tried and gave the same (2.04 s); the
  sweep is kept. The serial Debug `nitcad verify` total is taken at the close-out (step 14), against CV-Q19's 75 s.
- **Mutation:**
  - `energy_layer` (13 mutants: κ, G, the quartic's four factors, H's sign, the convergence threshold, the
    slowest-root choice, λ's factor, the scaling, a determinant term, the Poisson entry): 11 caught. The 2
    survivors are equivalent:
    - the scaling exponent (any scale gives the same roots);
    - the determinant term's sign: it multiplies the Poisson row's zero entry.
  - TVB-13 (10 mutants: the holes' term, the tail's 13, the rise's 2/3, the predicate, Canali's β, the drift
    criterion, the rise's denominator, the reference self-check, the Newton tolerance, the T_n criterion): 10
    caught. They are judged by any change in the printed comparison, failure and reason lines, because the full
    case is red.

**DT-Q37 (TVB-13's T_n reference).**
- (A) Evaluate the homogeneous reference at each judged node's own field (its two edges' mean) instead of V/L. The
  comparison is then the one 16.2 (f) describes. Measured as a diagnostic: 1.3e-14 … 1.5e-13, so all 4 would pass.
  The interior field's departure from V/L stays printed. This changes a reference input after a run (TVM2).
- (B) Record as found: 4 fails, pinned, with the cause printed.

**Recommendation: A.** V/L is wrong for this device: the contact layers take part of the voltage. 16.4 names T_n(E)
at the field, and the criterion is unchanged.

**16.43 Step 12 as built: thermal_slab, TVB-27 and TVB-28. DT-Q38 open.**
- **`verify/reference/thermal_slab`** (PD, base-free; κ = κ₃₀₀ (T/300)^−a, Unit 23 T11's form):
  - `kirchhoff` / `kirchhoff_inverse` in dd; `slab_temperature` (u linear in x).
  - `self_heated_resistor`: (κT′)′ = −J²/(qμN) − (3/2)(k/q) J T′ with n = N, T = T_amb at both ends and
    J = −qNV/∫dx/μ. Solved in dd: Dormand–Prince on (u, u′, ∫dx/μ) at fixed steps, Brent on u′(0), the secant
    on J. Second method: Kirchhoff-form finite differences on 2¹⁴ and 2¹⁵ cells (Newton, double), Richardson.
    After the first run its κ integral was moved from dd to a double closed form (no outcome changed: 51 s to
    3.2 s for TVB-27, the difference still 1.1e-11 K).
  - `seebeck_rise`: the closed form for Nc ∝ T^{3/2}, with Gauss–Kronrod as the second method (7.7e-34 eV);
    `peltier_coefficient` = −(kT/q)((5/2 + r) − ln(N/Nc(T))).
- **L0** (`thermal_slab_test`, 5 cases, 29 assertions): the Kirchhoff round trip (1e-29 K) and du/dT = κ; the
  slab's ends and sag (T(L/2) < 350 K); the constant-coefficient resistor against its closed form in dd (1e-12 K)
  and its steps-halving; the Seebeck second method and sign, the Peltier coefficient; the refusals. Two L0
  expectations of mine were wrong in the first draft (the sag's direction; the closed form evaluated in double
  lost 1e-9 K to cancellation) and were fixed before any NiTCAD run.
- **TVB-27 and TVB-28** (`tvb27.cpp`, Tier full, registered after TVB-19; 52 benchmarks). 1D ladder N = 512, 256,
  128, 64 with the (h, tol) factorial at N = 1024.
  - TVB-27: the slab (10 µm, 300/400 K ends, default silicon) and the self-heated resistor (10 µm, N_D = 1e17,
    4 V, default silicon). The resistor reference is solved once over every ladder's node set.
  - TVB-28 on the legacy silicon set (chosen before the run: Nc ∝ T^{3/2}, so S and Π are closed forms; the
    default Couderc set has no closed form): S as the open-circuit E_Fn rise on a 100 µm bar 300 → 400 K; Π as the
    odd part of the contact heat per current at ±1 mV; C5 on four states (the resistor at 4 V, the Seebeck bar,
    Peltier ±1 mV) with 11.13's Newton-term limitation printed.
  - **Construction change after a run error** (no outcome seen): Newton's 1e-12 tolerance with the factorial's
    tol/100 = 1e-14 did not converge at N = 512. tol_newton = 1e-10.
- **Frozen outcome, set before the first run:** every comparison passes.
- **TVB-27's run: 2 pass, 2 not_evaluated** (3.2 s Release). The frozen case fails, so `tvr4_full` fails, pending
  DT-Q38.
  - The slab is exact: errors 2.8e-16 at every level, so its order is 0 and "not in the asymptotic range".
    Accuracy passes.
  - The resistor's errors are flat at 1.66e-6, 1.55e-6, 1.14e-6 and 7.3e-7 relative (orders −0.09, −0.44,
    −0.63). Accuracy passes (≤ 1e-4); the order is not_evaluated.
  - **Cause, measured with temporary runs (not committed):** the floor is the reference's n = N assumption. The
    mobility gradient puts a space charge δn/N of about 4e-5 into NiTCAD's solution, scaling as 1/L². On the
    legacy set the floor halves (the Couderc T·M′/M Thomson term is the other half). On the legacy set at 100 µm
    the errors are 2.25e-6, 5.8e-7, 1.41e-7 and 2.9e-8, orders 1.93, 2.04 and 2.26: clean second order.
  - This is a construction error of the benchmark (n = N is this step's reference), not a finding.
  - An observation, not a finding: `tests/solve/electrothermal_test.cpp:652` says the 10 µm resistor's space charge
    is "some 1e-8" of N_D; at 4 V it is nearer 4e-5. Its tolerance holds either way.
- **TVB-28's run: 6 pass, 1 not_evaluated, 1 fail** (0.22 s Release).
  - S: accuracy 2.06e-9 passes. Its errors are 2.06e-9, 7.15e-9, 2.77e-8 and 1.11e-7, orders 1.79, 1.94 and 1.99
    (finest pair first). The finest pair is 0.01 below the window; not_evaluated.
  - **Fail: Π's order.** Π is exact in h (errors 4.721e-10 at every level, from Newton). That is the class of
    DT-Q36 (an h-exact quantity judged by an order).
  - Π's accuracy 4.7e-10 passes.
  - C5 passes on all four states: residual/bound 5.2e-5 (resistor), 9.7e-3 (Seebeck), 6.7e-6 and 6.6e-6
    (Peltier ±1 mV).
- **Time:** both are full tier; `tvr4_full` grows by about 3.4 s in Release.
- **Mutation:**
  - `thermal_slab` L0 (17 mutants: the transform's two 300s and κ's exponent sign, the Joule factor, the Thomson
    3/2, ∫dx/μ, J's relation sign, the mobility denominator, the dense evaluation, the halved-step shot, the slab's
    interpolation, Seebeck's A, its 3/2 and the integrand's Nc exponent, the Peltier factor, two refusals): 15
    caught. The survivors:
    - κ's exponent sign: the L0 resistor has constant κ (a = 0). Caught in TVB-27 (A1).
    - the V = 0 refusal: still refused, by the shooting's empty bracket.
  - TVB-27 (12 mutants, judged by the printed comparison, failure, reason and note lines since the frozen case is
    red): 9 caught (κ's sign, the FD Thomson term, the Richardson step, the slab's end, the voltage, the second
    method's gate, the relative scale, the accuracy criterion, the node positions). The survivors:
    - the FD Newton Jacobian's mobility term (equivalent: only the convergence rate changes);
    - u_iter = 0 (equivalent here: e(h, tol) = e(h, tol/100) exactly);
    - p_th 2 → 1 (not observable while both orders are not_evaluated).
  - TVB-28 (10 mutants): 9 caught (T1, the rise's offset, Π's T, Π's sign, the odd part, the GK gate, the material
    set, the thermopower exponent, C5's metric). The survivor: C5's outcome forced to pass, the same outcome as
    found (C5's bound is TVR-2's, mutated there).

**DT-Q38 (TVB-27 and TVB-28's non-passes).**
- (A) Record as found, as DT-Q36 A: the frozen cases pin the 8 passes, the 3 not_evaluated and Π's order fail with
  their reasons. The resistor's cause stays in the notes.
- (B) Construction change for the resistor only (TVM2): 100 µm on the legacy set, measured orders 1.93/2.04/2.26,
  so it would pass. The slab, S and Π stay as found (A).
- (C) B, plus judge the h-exact quantities (slab, Π) as exactness, as DT-Q36 B.

**Recommendation: B.** The resistor's floor is the benchmark's own reference assumption (n = N), like DT-Q37's V/L,
and the fix keeps the criterion. The slab and Π are DT-Q36's class and should follow its ruling.

**16.44 Step 13: TVB-32 on the TVR-4 devices. Scope and frozen expectations, written before code and before any run
(TVM2).**
- **Scope.** Each TVR-4 benchmark already built (TVB-05, 06, 10, 13, 15, 16, 17, 18, 19, 27, 28) exposes a
  `tvbXX_consistency()` that runs the checks on its own device (its own builder and models) at the coarsest level
  of its ladder.
  - The operating points are the benchmark's own, plus the zero-bias point where the benchmark has none.
  - TVB-20 and 21 join with step 9.
- **Checks:**
  - C1–C5 and C12 at every DC point (a small-signal run's DC points included);
  - C1 and C3 on the equilibrium solve, for devices with a semiconductor region;
  - C7 and C8 at every small-signal point;
  - C6 on every transient run (TVB-10's deep depletion, TVB-19's two cases);
  - C9 on TVB-18's diode at 0.40 V, from a DC stencil of h = 1 mV, with f = 0 among the frequencies;
  - C10 and C11 are not run: no TVR-4 benchmark builds a limit series or a symmetry map. Listed in the notes.
- **Tier.** TVB-32's fast tier is unchanged (the three decks). The full tier adds these entries. The frozen case runs
  in `tvr4_full` as `[.tvr4_full]`.
- **Frozen expectations** (11.13's rule, unchanged):
  - C1, C2, C7, C8: pass where every bias is exactly 0 and the thermal contacts are at one temperature, else NOT
    APPLICABLE naming the predicate.
  - C1 and C3 on the equilibrium: pass. C3, C4 and C12: pass. C6 and C9: pass.
  - C5: NOT APPLICABLE ("electrothermal on") with the model off; pass with it on (TVB-27, 28), with 11.13's
    Newton-term limitation printed.
  - A failed or NOT EVALUATED part, and any other outcome, is a finding to the owner (DT-Q11: not fixed in TVR-4).
- **Prediction, not an expectation** (written before the run): TVB-15's zero-bias point carries the TVQ3 hook's
  manufactured generation. C1's header says "NiTCAD has no optical or prescribed generation, so that predicate
  holds by construction", and the hook makes that untrue. C1 and C2 at that point are therefore expected to fail
  (the generated carriers leave through both contacts). If they do, that is a finding against C1's predicate list
  (11.6), not against the solver.

**16.45 Step 13 as built: TVB-32 on the TVR-4 devices. One finding; DT-Q39 open.**
- **Code.**
  - `tvb32.hpp` exposes TVB-32's collector (`Tvb32Collector`: dc, ac, transient, equilibrium) and one
    `tvbXX_consistency()` per benchmark, each in its benchmark's file on its own builder and models.
    `tvb32_tvr4_entries()` runs them in order.
  - `run_tvb32(Tier)`: fast is unchanged (the three decks); full and nightly add these entries.
  - To share their builders, four benchmarks had device or case construction moved into functions, with no
    change in what they compute: TVB-05 `coax`, TVB-06 `stack`, TVB-10 `dd_transient` and `flat_band`, TVB-15
    `generated_point`, TVB-19 `relaxation_case`, `pulse_case`, `case_options` and `solve_case`.
  - **Checked as identical:** a baseline build of HEAD (be65856) in a separate worktree. TVB-05, 06, 10 and 15
    (`-s`), TVB-19's nightly case and TVB-32's deck test (7037 lines) print the same output, sorted, line for
    line. Catch2's case order differs between the builds.
- **Devices and points** (the coarsest level of each ladder):
  - TVB-05: coax M = 32. TVB-06: 1D N = 32. Both at (0, 0) and (0, 1) V, 1 MHz.
  - TVB-10: N = 64. HF (equilibrium Poisson) at 0, V_FB + 0.5 and V_FB + 2 V, three frequencies; equilibrium;
    the deep-depletion ramp.
  - TVB-13: the drift bar N = 64 at 0 and four fields; the energy-balance bar N = 128 from 0 to 1e4 V/cm by half
    decades; equilibria.
  - TVB-15: N = 64, zero bias without and with the E part's generation; equilibrium.
  - TVB-16: the three configurations at N = 128, forward from 0 and reverse; equilibria.
  - TVB-17: N = 64, AC at 0, −0.5, −1 and −2 V.
  - TVB-18: N = 128, AC at 0, 0.35, 0.40 and 0.45 V at 0, 1e4 and 1e6 Hz; C9 at 0.40 V (h = 1 mV).
  - TVB-19: both cases at N = 257 and 8 dt₀ (relaxation BDF2 and BE, the pulse BDF2), with each bar's zero-bias
    point and equilibrium.
  - TVB-27: the slab and the resistor (0, 4 V) at N = 64. TVB-28 (legacy set): the Seebeck bar, the Peltier bar
    (0, ±1 mV) and the resistor (0, 4 V), all at N = 64; equilibria.
- **Run: 457 entries in 0.53 s (Release): 261 pass, 175 not applicable, 19 NOT EVALUATED, 2 fail.** 52 assertions
  of the frozen case fail, so `tvr4_full` stays red pending DT-Q39.
  - Every check where the frozen rule applies is as frozen. The worst ratio is 0.073 of the bound (C8, TVB-10 HF).
- **Departures from 16.44's expectations, classified:**
  - **The finding: C1 and C2 fail at TVB-15's generated zero-bias point.** This was predicted in 16.44.
    - Edge currents 8.4e-9 against 5.4e-12; quasi-Fermi split 0.22 eV against 1.3e-14 eV; SRH 4.3e-8 against
      5.4e-12.
    - The TVQ3 hook (`models.manufactured`) is a prescribed generation, so C1's "no optical or prescribed
      generation, holds by construction" is false. `DeviceData` and `BiasPoint` do not record the hook, so C1 and
      C2 cannot see it.
    - The checks are wrong to apply here, not the solver: the currents are the generated carriers leaving.
    - Not fixed (DT-Q11).
  - **Departures that are the checks' own contracts (11.6), not findings.** The frozen rule was 11.13's, written
    for decks with every mechanism on and no insulator-only device:
    - C2 is NOT EVALUATED, "no rates stored", at the zero-bias points with G-R off (9 entries).
    - On insulator-only devices (TVB-05, 06), C1 and C2 are NOT EVALUATED ("no nonzero terminal resolution") at
      zero bias, and C12 ("no semiconductor nodes") at every point (8). DT-Q2 A made these devices after C1–C12
      were built.
    - At TVB-10's zero-bias HF point (equilibrium Poisson: no currents), C1 is NOT EVALUATED, "no nonzero terminal
      resolution" (1).
    - C9 at TVB-18 is NOT EVALUATED: its conductance part passes (3.2e-4 of the bound), and its capacitance part
      needs a gate or electrode row (1).
    - C6 is NOT APPLICABLE on TVB-19's three transients ("transient with a gate or electrode"); TVB-10's ramp
      passes (3).
- **Mutation** (12 mutants, judged by any change in the printed tally, failure and reason lines, since the frozen
  case is red): 11 caught. These were the thermal-temperature predicate, removing C12, C8, the equilibrium C3,
  TVB-18's entries or C6, C9's step, TVB-15's generation, the slab's temperatures, TVB-13's energy model and the
  zero-bias flag. The survivor, C6's device data from the HF models, is equivalent: C6 reads none of the models
  that differ.

**DT-Q39 (TVB-32 on the TVR-4 devices).**
- (A) Pin as found:
  - the frozen case's expectations add the five contract classes above, each with its reason;
  - the two C1/C2 fails at TVB-15's generated point become a known finding;
  - the finding's fix (a "prescribed generation" predicate from `models.manufactured` into `DeviceData`, making
    C1 and C2 NOT APPLICABLE there) is deferred past TVR-4 (DT-Q11).
- (B) As A, but authorize the predicate fix now as a TVR-4 change (TVM2's change to a check's inputs; its own
  tests, bit-identical elsewhere).

**Recommendation: A.** DT-Q11 keeps findings out of TVR-4. The contract classes are what C2, C6, C9 and C12 state.

**16.46 Step 9: pao_sah, TVB-20 and TVB-21. The reference and the construction, recorded before any NiTCAD run of
the benchmark (TVM2).**
- **`verify/reference/pao_sah`** (PB through the two reproductions of 16.17; base-free):
  - The field E(ψ, V) in both forms of 16.17, with exact bulk densities (`neutral_bulk`).
  - ψ_s(V) by Brent in dd from V_G − V_FB = ψ_s + εE(ψ_s, V)/C_ox.
  - Q_I(V) by adaptive Gauss–Kronrod in dd from c ψ_B to ψ_s; Q = μ ∫₀^{V_D} (−Q_I) dV by Gauss–Kronrod.
  - **The primary reference is form 1 with c = 1** (Wong's reproduction, the one that states its limits). The
    other form and c = 1/2 enter only as ε_cutoff (16.18 (d)).
  - At c ψ_B > 0 the inner integrand is finite (16.18 (n)).
  - The second method (16.18 (d)): g_d against a fourth-order difference of Q in V_D at the dyadic step 2⁻¹⁷. It
    is run at L0 and at TVB-20's run time; above 1e-10 it is a run error.
  - Brews' charge sheet (16.18 (e)), with ψ_s at both ends from the same gate equation instead of Brews' χ-modified
    eq 5. A diagnostic only; its closed form is checked against Gauss–Kronrod.
  - Quadrature tolerances: 1e-24 (inner) and 1e-22 (outer), relative. The quadrature estimate enters U_ref.
- **L0** (`pao_sah_test`, 5 cases, 62 assertions):
  - at V = 0, ψ_s against TVB-09's Kingston–Neustadter reference (separate code), 1e-24 relative, both forms;
  - the bulk ambiguity in strong inversion: the forms within 1e-10 of Q, the cutoffs within 1e-8;
  - the identity at weak and strong inversion, 1e-12;
  - the linear-region limit, 1e-3;
  - the subthreshold slope against the depletion formula, 5 %;
  - the charge sheet within 5 % in strong inversion, and its second method, 1e-24;
  - the refusals.
- **Corrections to my own L0 expectations before any NiTCAD run:**
  - The identity needed a fourth-order stencil: the central difference's (h/V_T)²/6 is 2.5e-6.
  - It needed a dyadic step, because the double rounding of V_D ± h gave a 2.7e-12 floor.
  - The cutoff bound was 1e-10 by my guess; the estimate n_i V_T / (E(ψ_B) N_inv) ≈ 2e-9 (measured 5.6e-10)
    makes it 1e-8.
- **L0 time:** the ambiguity and identity cases take 1.5 s in Release, against CV-Q19's Debug budget. They are
  hidden (`[.pao_sah_release]`) and run in the Release-only ctest `verify.pao_sah_release`, as `diode_profile`
  does. The identity also runs at TVB-20's run time.
- **Construction** (`tvb20.cpp`; 16.2 (i), 16.18 (d)):
  - **Device.** A uniform p-substrate, N_A = 1e17. Abrupt n+ source and drain boxes, 1e19, 0.5 µm wide and
    0.1 µm deep. Depth 1 µm, with the body contact on the bottom.
  - **Gate.** A lumped gate over the channel: t_ox = 10 nm, ε_ox 3.9, metal 4.5 eV, no fixed charge.
  - **Physics.** 300 K, default silicon. Constant mobility: doping, field and surface mobility off, so μ_n is
    μ_max. Boltzmann, no BGN, SRH, Auger and radiative off. Weighted terminal current (Unit 32).
  - **Lengths.** L = L0 × {1, 2, 4, 8} with L0 = 2.5 µm (W_dmax ≈ 0.11 µm, so L0 ≥ 20 W_dmax).
  - **Mesh, x.** Spacings 2 nm at the junction edges, growing ×1.2 up to 50 nm in the source and drain and up to
    L/40 in the channel.
  - **Mesh, y.** 0.2 nm at the surface, growing ×1.15 up to 40 nm.
  - **Factorial.** At 8 L0: tol 1e-12 and 1e-14, and h/2 (every spacing halved, the ratios square-rooted).
  - **U_disc** is the two-grid GCI (F_s 3, p 2) of h against h/2. U_iter = |Q(tol) − Q(tol/100)|.
  - **Bias points.** Q = |I_D| L.
    - V_D = 50 mV at V_G = 0.45, 0.6, 0.8, 1.2, 1.6 and 2.0 V.
    - V_D = 0.5 V at V_G = 0.6, 0.8, 1.2, 1.6 and 2.0 V.
    - Each V_D is swept in V_G from equilibrium.
    - A point is NOT APPLICABLE where the reference's ψ_s(V_S) < ψ_B + 3V_T.
  - **TVB-21.** S at V_G = 0.4, 0.5 and 0.6 V, V_D = 50 mV, from the stencil V_G ± 20 mV. The same stencil is
    used on NiTCAD's ladder and on the reference. The depletion formula is printed (DT-Q6).
  - **Judging.** The order in 1/L comes from the ladder's three differences (finest pair), with the spread of its
    two values as the asymptotic check. Q_∞ is Richardson's at the finest pair, and U_L its three-grid GCI.
  - **The accuracy needs the order.** It is NOT EVALUATED unless the order is a pass (Q_∞ is not established
    otherwise), or when ε_approx > 1e-2 Q_PS.
  - V_T and g_m are named in §6's catalogue row but not in 16.19's frozen criteria, so they are not judged.
- **Tier:** nightly (`verify_full.tvb20`, NITCAD_VERIFY_FULL_TIER, as TVB-19). Registered in full, like TVB-19.
- **Frozen outcome, set before the first run:** every comparison passes, except those at points below ψ_B + 3V_T,
  which are NOT APPLICABLE by that predicate. 22 comparisons for TVB-20, 6 for TVB-21.
- TVB-32 adds TVB-20's device at L0 (the V_D = 50 mV sweep from zero bias, and the equilibrium) under 16.44's
  rule.

**16.47 Step 9 as built: TVB-20 and TVB-21. DT-Q40 to DT-Q43, ruled A ("for decision choose recommended",
2026-11-11).**
- **First attempts: run errors, no comparison seen.**
  - With 16.46's tol 1e-12 (and an interim 1e-10) Newton did not converge at L = 2 L0, V_D = 0.5 V, V_G = 0.8 V, with
    the bias step cut to 2⁻³⁰. The last updates wander between 2e-10 and 8e-10 without decreasing.
  - **Cause, measured with a temporary print in `DriftDiffusion::update_size` (not committed):** the largest
    correction is always a hole density at a node of the inverted or depleted surface. There p is 7e-22 to 4e-19 of
    the largest density (1e-16 of the local electrons) and its correction is 1e-10 of itself: the linear solve's
    rounding. Unit 50's `minority_floor` does not help (ε n is below `update_size`'s own 1e-20 floor there).
  - These holes carry no charge or current that enters I_D. The solver is not changed.
- **DT-Q40 (the Newton tolerance). A, applied:** tol 1e-6, so the factorial's tol/100 is 1e-8, the solver's default,
  12 times above the worst stall seen. U_iter came out at 1e-13 of Q or below at every point, so the choice enters
  no outcome. `solve_level` now prints the stopped point's last six updates in a run error.
- **First complete run (L0 … 8 L0, 8.8 min Release): TVB-20 16 pass, 6 NOT EVALUATED; TVB-21 1 pass, 4 NOT
  EVALUATED, 1 fail.**
  - TVB-20: the 11 orders pass (0.99 to 1.08). Five limits pass. Six are NOT EVALUATED: U_L is 1.5 to 2 % of Q_PS
    at V_G ≤ 0.8 V (and at 1.2 V, V_D = 0.5 V), above 16.19's 1e-2.
  - TVB-21: S at 0.6 V missed its bound by 1.19 (1.67e-5 V/dec against 1.40e-5, 2e-4 of S); at 0.4 and 0.5 V the
    spread of the two orders was 0.38 and 0.11.
- **DT-Q41 (U_disc's form). A, applied.** `tvb20.cpp` computed F_s |f(h/2) − f(h)| / (r^p − 1), the estimate of the
  h/2 value. The judged values are on h. 16.36's form (`registry.cpp`), F_s |f(h/2) − f(h)| r^p / (r^p − 1), is the
  one 16.46 names ("the two-grid GCI of h"). A coding error of mine, found after the run from the S fail; four times
  larger. It enters the bounds only, never ε_approx.
- **DT-Q42 (the six NOT EVALUATED limits). A, applied: the ladder is L0 × {4, 8, 16, 32}** (10 to 80 µm), the
  factorial at 32 L0. The criterion is unchanged. (B was to pin the six as found.)
- **The ruled run (12 to 16 min Release, two runs): TVB-20 22 pass. TVB-21 6 NOT EVALUATED.**
  - TVB-20: orders 1.006 to 1.085, spreads within 0.1. U_L is 0.06 to 0.52 % of Q_PS. |Q_∞ − Q_PS| / bound is 0.08
    to 0.27.
    - Q_∞ against Q_PS: −0.37 to −0.68 % at V_G ≤ 0.8 V, where U_disc is 1.8 to 2.1 %; within ±0.08 % at V_G ≥ 1.2 V,
      where U_disc is 0.03 to 0.23 %. The weak-inversion bound is the mesh's, not the ladder's.
    - No point is below ψ_B + 3 V_T (0.495 V; ψ_s is 0.549 V already at V_G = 0.4 V).
    - Brews' charge sheet at V_G 1.2 V, V_D 50 mV: 2.5 % below Pao–Sah (diagnostic).
  - TVB-21: S does not depend on L measurably at these lengths. Its ladder's differences are 3e-7 to 5e-5 V/dec,
    below U_disc (4.4e-5 to 5.6e-5), so the orders are −1.58, −0.23 and 0.30 with spreads 0.6 to 1.2: "not in the
    asymptotic range", and the limit "Q_inf not established".
    - Printed, not judged: S at 32 L0 against the reference is +0.84 % at 0.4 V (ε_cutoff 0.75 %), +0.069 % at 0.5 V
      and +0.021 % at 0.6 V; on h/2, +0.82 %, +0.053 % and +0.007 %.
    - The depletion formula is 2.5 to 2.7 % below the reference's S (DT-Q6's diagnostic).
- **DT-Q43 (TVB-21's outcome). A: pinned as found**, as DT-Q20 A and DT-Q35 A. (B, a criterion change, was to judge
  S at the longest L with ε_approx from the last ladder difference when the differences are below U_disc. It is
  DT-Q36's class, an L-exact quantity judged by an order, and should follow that ruling if B is ever taken.)
- **Frozen case** (`[tvb20_full]`, `verify_full.tvb20`): TVB-20 22 pass; TVB-21 6 NOT EVALUATED with the two
  reasons. 16.46's "every comparison passes" is replaced.
- **TVB-32 on TVB-20's device** (L0, zero bias and the V_D = 50 mV sweep, the equilibrium): 44 entries, 501 in all.
  C1 2 pass and 6 not applicable, C3 8 pass, C4 7 pass (worst 0.0097 of the bound), C12 7 pass, C5 7 not
  applicable, C2 6 not applicable and 1 NOT EVALUATED ("no rates stored" at zero bias with G-R off: 16.45's first
  contract class, now 10 entries). No finding. The frozen case stays red pending DT-Q39.
- **Mutation.**
  - `pao_sah` L0 (19 mutants: the field's scale, ψ_B, the hole term's sign, each form's electron term, n's V, the
    gate equation's sign, the cutoff, Q_I's sign, μ, the outer integral's lower limit, the stencil's 12, the
    identity's sign, the charge sheet's V_T, its 2/3 and its diffusion sign, three refusals): 16 caught at first
    (the hole term's sign by a 120 s limit: the root search does not end).
    - Two survivors were gaps and are closed: the cutoff ignored (now: c = 1/2 must change the current) and the
      charge sheet's ψ_s ≤ V_T refusal (now tested). `pao_sah_test` is 65 assertions.
    - One is equivalent: without the V_G ≤ V_FB refusal the point is still refused, by Brent's empty bracket.
  - TVB-20 and TVB-21's own code: not run. Each mutant costs the 12 to 16 min nightly run.
- **Time:** nightly only; `verify_full.tvb20` 12 to 16 min Release against its 3600 s limit. Nothing is added to the
  fast tier or to `tvr4_full` except TVB-32's 44 entries (the L0 sweep).

**16.48 DT-Q35 to DT-Q39 ruled as recommended ("for decision choose recommended", 2026-11-11) and applied.
`tvr4_full` is green.**
- **DT-Q35 A (TVB-10): pinned as found.** 13 pass; C_HF at V_G − V_FB = 2 V is NOT EVALUATED ("U_num not
  estimated"). The five C_HF passes' bounds are pinned to 1e-6. Nothing else changed.
- **DT-Q36 A (TVB-19, nightly): pinned as found.** 10 pass, 3 NOT EVALUATED ("not in the asymptotic range": the
  decay rate's and the width's BDF2 orders in dt, the width's order in h), 3 fail (the orders in h of the decay
  rate, the pulse's area and its centre, which have no h error). The run of 2026-11-11 reproduced 16.41's outcome
  line for line. Its criteria carry no computed bound, so none is pinned.
- **DT-Q37 A (TVB-13): the T_n reference is taken at each judged node's own field** (its two edges' mean), as
  16.2 (f) states, not at V/L.
  - T_n − T_L: 3.7e-14 and 1.3e-14 (N = 128), 1.5e-13 and 6.8e-14 (N = 512), against 2.3e-10 and 3.7e-9. All 8
    comparisons pass; the case is `require_all_pass` again.
  - Printed, not judged: the rise against the homogeneous rise at V/L (the former comparison) and the interior
    field's departure from V/L.
- **DT-Q38 B (TVB-27): the self-heated resistor is 100 µm on the legacy set.** The slab, S and Π stay as found.
  - The resistor's errors are 2.25e-6, 5.83e-7, 1.41e-7 and 2.93e-8 (N = 64 … 512): orders 1.93, 2.04 and 2.26.
    Accuracy passes (2.9e-8 against 1e-4).
  - **Its order is still NOT EVALUATED**: p̂ is 2.26, inside the window, but the successive orders spread 0.225,
    above 0.1. 16.43's "so it would pass" overlooked the spread rule; the diagnostic's numbers were these.
  - TVB-27 as pinned: both accuracies pass, both orders NOT EVALUATED. TVB-28 as pinned: 6 pass, S's order NOT
    EVALUATED, Π's order fail (DT-Q36's class).
  - TVB-32's TVB-27 resistor is the same device (100 µm, legacy set): C1, C3, C4, C5 and C12 pass.
  - TVB-28's C5 state "self-heated resistor, 4 V" is unchanged (10 µm, legacy set).
- **DT-Q39 A (TVB-32): pinned as found, 501 entries.** The frozen rule of 16.44 stands. Beside it the case pins,
  by count and reason:
  - C2 NOT EVALUATED, "no rates stored" (10);
  - C1 and C2 NOT EVALUATED, "no nonzero terminal resolution" (5); C12 NOT EVALUATED, "no semiconductor nodes" (4);
  - C9 NOT EVALUATED at TVB-18, no gate or electrode row (1);
  - C6 NOT APPLICABLE, "transient with a gate or electrode" (3);
  - the known finding: C1 and C2 fail at TVB-15's generated point (2). Its fix stays deferred past TVR-4 (DT-Q11).
  Any other departure fails the case.
- **Tests.** `device_theory_test` gains `require_as_found` (quantity, outcome, a part of the reason, and the bound
  to 1e-6 where the criterion states one). `[.tvr4_full]`: 11 cases, 1041 assertions, all pass (Release).
- **Not done:** mutation runs of these changes (the T_n reference's field, the resistor's set and length, the
  pins) and of TVB-20/21's own code (16.47); the serial Debug timing (CV-Q19's 75 s), which these full-tier
  changes do not enter.
- **Open in TVR-4:** step 14, the close-out. No DT question is open.

**16.49 DT-Q44 and DT-Q45, ruled A (owner "go ahead and choose recommended", 2026-11-11).**
- **DT-Q44 (TVB-21 judged nothing after 16.47).**
  - (A) Where the last ladder difference of S is within U_disc, S has no measurable dependence on L: its order in
    1/L is NOT APPLICABLE by that predicate, and S at the longest L is judged directly,
    |S(32 L0) − S_PS| ≤ U_num + U_ref + ε_approx with ε_approx = |last ladder difference| + ε_cutoff. Elsewhere
    16.19's rule stands. This is a criterion change after a run (TVM2), DT-Q36 B's class.
  - (B) Leave the six NOT EVALUATED (16.47's DT-Q43 A).
  - **A, applied.** Without it the swing benchmark verifies nothing.
  - **Expectation, set before the run from 16.47's printed ladders:** all three orders NOT APPLICABLE; S passes at
    0.97, 0.62 and 0.36 of its bound (0.4, 0.5 and 0.6 V).
  - **Run: as expected.** 0.964, 0.617 and 0.360. The last differences are 4.6e-5, 4.0e-6 and 3.5e-7 V/dec against
    U_disc 5.6e-5, 4.9e-5 and 4.4e-5.
  - The pass at 0.4 V is marginal. Its bound is 86 % ε_cutoff: at ψ_s − ψ_B = 5 V_T the reference's cutoff
    ambiguity is 0.75 % of S, and NiTCAD's S lies inside it. At 0.5 and 0.6 V the bound is mostly U_disc.
  - TVB-20 is unchanged: 22 pass. Frozen case: TVB-20 22 pass; TVB-21 3 NOT APPLICABLE, 3 pass.
- **DT-Q45 (CV-Q19's 75 s after TVR-4).**
  - Measured, serial Debug `ctest -R verify -LE full`: 71/71 pass in 101.3 and 94.9 s, **under load** (a game on
    the machine). The 66 tests of 16.38's idle run (74.5 s) ran 22 to 30 % slower than there.
  - Idle estimate from the per-test times: 71.05 s for the 62 tests unchanged since 16.38, plus about 6.2 s for
    the changed and new ones, about 77.2 s. Over the budget by about 2.2 s. The growth is TVB-13's fast case
    (about 1.4 s; 16.42) and the L0 tests of steps 8 to 12 and 9.
  - (A) Move the Debug-heavy L0 cases to a Release-only ctest, as 16.31, 16.38 and 16.46 did.
  - (B) Raise CV-Q19 to 80 s.
  - **A, applied.** Eight cases carry `[.tvr4_reference_release]` and run in the Release-only ctest
    `verify.tvr4_reference_release`: `pao_sah` limits; `admittance` at ω = 0, the Riccati second method and C_SCR;
    `diode` SNS; `minority_diffusion` erfc and collocation; `thermal_slab`'s constant-coefficient resistor.
    Nothing is removed from Release or CI.
  - After the move: 71/71 pass in 97.0 and 96.7 s, under the same load (the unchanged tests 30 % slower than
    idle). **Idle estimate 74.2 s.** This is an estimate; no idle run was possible today.

**16.50 Step 14: TVR-4 close-out (2026-11-11).**
- **Outcomes** (ARCHITECTURE.md 10 carries the table): TVB-06, 13, 15 and 20 pass in full. TVB-21 3 pass and 3 NOT
  APPLICABLE. TVB-05, 10, 16, 17, 18 and 27 pass where judged and are NOT EVALUATED elsewhere (TVB-16 also 3
  finding, (x)). TVB-19 and TVB-28 carry 3 and 1 order fails of h-exact quantities. TVB-32 on these devices: 501
  entries, 285 pass, 194 not applicable, 20 NOT EVALUATED, 2 fail (the finding of 16.45).
- **Provenance:** PD 9 references (TVB-05, 06, 10, 13, 15, 17, 19, 27, 28), PA 2 (TVB-16, 18), PB 1 (`pao_sah`:
  TVB-20, 21; Pao and Sah's paper was not read, 16.17).
- **Mutation, for 16.47 to 16.49's changes.**
  - TVB-20 and TVB-21's own code (12 mutants, each a nightly run): 11 caught. These were U_disc's r^p, ε_cutoff
    out of the flat bound, the flat predicate, the ladder's first length, Q without L, the 1e-2 limit, the order
    window, the reference's μ by 1 % (6 comparisons fail), ε_cutoff's first term, the reference's S stencil and
    the spread limit.
    - The survivor: S judged at Q_∞ instead of at the longest L. Same outcomes: the two differ by less than the
      bound's slack.
  - TVB-13 (3): the node field by 0.1 % and the reference at V/L are caught. The outcome forced to pass is the
    outcome as found.
  - TVB-27 (3): the resistor back at 10 µm survived the first pins, since both lengths pass ≤ 1e-4. The case now
    pins the finest error below 1e-7 (2.9e-8 as found; about 1e-6 at 10 µm), which catches it. The default set's
    κ and μ in the reference are equivalent: the default set copies the legacy's mobility and conductivity.
  - Not mutated: the pins of 16.48 (test code).
- **Suites (2026-11-11, after every change above).**
  - Release, `ctest -LE full -j 4`: 206/206 pass in 807 s. `verify_benchmark.tvr4_full` 77.7 s under load, against
    16.7's 120 s; `verify.tvr4_reference_release` 0.49 s; `verify.pao_sah_release` 1.54 s.
  - Debug, serial `ctest -R verify -LE full`: 71/71 pass (16.49's times).
  - Nightly cases run today in Release: `verify_full.tvb20` passes (38 assertions, 13.6 min) and
    `verify_full.tvb19` passes (25 assertions, about 9 min). Not rerun: `verify_full.tvb04`, `tvb05`, `tvb07`,
    `tvb09` and `tvb14`, whose code this branch's last three commits do not touch.
- **Time.** Fast tier: TVB-06 0.09 s and TVB-13 1.65 s in Debug under load (about 1.3 s idle), within 16.7's 2 s
  only at idle. Serial Debug verify: 16.49. Nightly: TVB-19 about 9 min, TVB-20/21 12 to 16 min.
- **Not done, for the owner:**
  - an idle serial Debug run to confirm the 74.2 s estimate (CV-Q19);
  - findings deferred past TVR-4 (DT-Q11): the prescribed-generation predicate for C1 and C2 (16.45); finding (x);
  - DT-Q36 B's redesign (judging h-exact quantities as exactness), which would turn TVB-19's and TVB-28's four
    fails and TVB-27's slab into judged comparisons;
  - TVB-16, 17 and 18's NOT EVALUATED comparisons (U_disc outside the asymptotic range) stay as found.
- **Work location.** Steps 9 and 14 were done on `verify/device-theory-step9`, a worktree branch off
  `verify/device-theory` at 7688cd6. Nothing is pushed. TVR-4 is complete on that branch.
