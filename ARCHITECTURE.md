# NiTCAD — Architecture

**PROPOSAL, accepted as the current proposal. The architectural questions needed for Unit 1 are decided (section 14). The owner named Unit 1 (build scaffold) on 2026-10-02.**

The owner accepted this architecture as the current proposal and recorded decisions D1–D6, D8 and R1–R5.
Section 14 classifies each as decided, provisionally decided, deferred or blocking verification. On 2026-10-02 the owner
accepted the advisor recommendations in `DECISIONS.md` for Q1–Q5, Q10 and Q11 (14.1). It is still a proposal: no
unit is complete, and the document must be rewritten to describe what actually exists as units land.

Naming: **NiTCAD** is this project. **NT-SemTCAD** is the legacy reference repository
(`C:\Users\disha\OneDrive\Desktop\NT-SemTCAD-claude-zealous-ritchie-kuuwr9\NT-SemTCAD`), used only as a
source of science, algorithms, tests and benchmarks. Legacy paths below are relative to its `pytcad/`
directory. Legacy facts in this document are tagged **[verified]** (checked against the legacy files named),
**[derived]** (computed here from verified inputs) or **[unverified]** (not checked; do not rely on it).
Section 12 lists them.

---

## 1. Fixed constraints (already decided by the owner)

| Item | Decision |
|---|---|
| Platform | Windows only |
| Languages | C++23 (primary); C and Fortran only where technically justified |
| C++ standard | C++23, enforced by the build system, not by compiler defaults |
| Desktop layer | Win32 + Direct3D 12 + Direct2D. No Qt, VTK, Python GUI or other GUI framework |
| Runtime | No Python |
| Tests | Written in C/C++ |
| Git | One branch per architectural unit; no auto-merge; `main` stays clean |
| Scope | Feature freeze: only explicitly requested things are built |
| Physics | Preserved: equations, units (cm, cm⁻³, V, A·cm⁻², s), constants, sign conventions |

## 2. Design principles

1. One responsibility per file; explicit dependencies; no global state; no circular dependencies.
2. Dependencies point one way (section 3). A lower layer never includes a higher one.
3. **Physics, discretization, solving and device description are separate.** The legacy
   `device.py`, `device2d.py` and `device3d.py` are large per-dimension classes (80 K, 120 K and 123 K) that
   hold geometry, scaling, physics, assembly and solve together. [verified: file sizes]
4. **Dimension is a parameter, not a code path.** The mesh is one node/edge/control-volume graph (box
   method) valid for 1D, 2D and 3D, so no unit is written for one dimension and later rewritten.
   Dimensional-reduction identities (uniform 3D → 2D → 1D) are test gates.
5. The solver core knows nothing about the UI. The UI knows nothing about physics. They meet only
   through plain data (device description in, results out).
6. Architecture changes do not change physics. Tolerances, not bit-identity, define agreement (section 6.8).

## 3. Layers and dependency diagram

Each layer may depend only on layers drawn below it. Layers on the same row do not depend on each other.

```text
 L8   app         Win32 application: windows, input, run control UI
 L7   render      Direct3D 12 / Direct2D drawing of meshes, fields, plots
 L6   analysis    derived observables from results   (deferred; nothing requested yet)
 L5   solve       Newton, damping, continuation, sweeps; produces results; cancellation + progress
 L4   assemble    scaling; discretization (SG flux and extensions); contacts/BC handling;
                  residual + Jacobian; nonlocal models (deferred)
 L3   device | results      device = problem description (regions, doping, materials, contacts)
                            results = plain result data (fields, terminal quantities, run record)
 L2   physics     materials and local models; value + exact partial derivatives; no mesh, no matrices
 L1   linalg | mesh         linalg = sparse matrix + direct solve   mesh = node/edge/control-volume graph
 L0   base        constants (CODATA 2018), units, error types
```

Allowed dependencies (everything not listed is forbidden):

| Layer | May depend on |
|---|---|
| base | nothing |
| linalg | base |
| mesh | base |
| physics | base |
| device | base, mesh, physics (material parameters only) |
| results | base, mesh |
| assemble | base, linalg, mesh, physics, device |
| solve | base, linalg, assemble, device, results |
| analysis | base, mesh, results |
| render | base, mesh, results, analysis |
| app | all of the above except physics and assemble |

`results` sits below `solve` so that `analysis`, `render` and `app` can read results without including any
solver header.

The layer names above are decided (D5). The root namespace is `NiTCAD`. Nesting each layer as
`NiTCAD::<layer>` is the working assumption and is provisional (Q1).

## 4. Layer responsibilities

| Layer | Owns | Does NOT own |
|---|---|---|
| base | physical constants, unit convention, error types | physics |
| linalg | sparse matrix storage, backend-neutral solver interface (6.10), Eigen SparseLU backend (D3), residual check | meshes, Newton |
| mesh | node/edge/control-volume graph, per-edge length and coupling area, per-node volume; tensor-grid constructors for 1D/2D/3D; no tensor-grid indexing in its public interface (6.11) | doping, materials, linear algebra |
| physics | material parameters, local models (mobility, recombination, …) | meshes, matrices, bias, contacts |
| device | regions, doping, material assignment, contact definitions (data only) | solving, assembling |
| results | field and terminal-quantity containers, run record, convergence history | solving |
| assemble | scaling, Scharfetter–Gummel flux (and its extensions), boundary/contact residuals, F(x) and J(x), gathering node/edge inputs for models | Newton iteration, UI |
| solve | Newton, damping, continuation, sweeps, cancellation, progress | assembly, rendering |
| analysis | quantities derived from results | solving |
| render / app | pixels, windows, user interaction | numerics |

The legacy placed the SG flux inside discretization together with physics hooks. Here it is an assembly
concern, because its extensions couple to physics (legacy handles heterojunction band offsets and
thermionic emission through per-edge arrays in the `Inputs` struct, `core/include/tcad/device1d/inputs.hpp`
[verified]).

## 5. Physics model contract

1. A model is a pure function of explicit inputs. It holds no mesh, no matrix and no global state.
2. A model returns its **value and the exact partial derivatives** with respect to each input the Newton
   system differentiates by. The legacy `RecombinationResult` in `core/include/tcad/physics/kernels.hpp`
   follows this pattern: `recombination_boltzmann()` returns R, dR/dn and dR/dp [verified: struct and function opened]. Every model's derivatives are
   validated against finite differences.
3. **Node-local models** (SRH, Auger, bulk Caughey–Thomas, band-gap narrowing) take values at one node.
4. **Edge-local models** (field-dependent mobility, impact-ionization driving field) take values gathered
   by the assembler for one edge, or for one node from its incident edges. The legacy local impact
   ionization does exactly this: per node it averages the incident edges per axis
   (`pytcad/ii_grid.py` docstring [verified]). The assembler gathers inputs and scatters the chain rule; the
   model never sees the mesh.
5. **Nonlocal models** (path-integrated band-to-band tunnelling, nonlocal ionization) need paths across the
   mesh, so they live in `assemble`, not `physics`. They are deferred.
6. The first physics unit contains only what the first diode needs (section 11).
7. **Carrier statistics (R4):** Boltzmann statistics initially; Fermi–Dirac and others are deferred. Statistics
   stay isolated inside `physics` so that adding Fermi–Dirac later is an addition, not a signature break for callers.
   Evidence from the legacy: `recombination_fd()` additionally takes the equilibrium product `np_eq` and its
   partials `dnpq_dn`, `dnpq_dp`, which `recombination_boltzmann()` does not
   (`core/include/tcad/physics/kernels.hpp` [verified]). The Unit 5 signatures are to be designed with that in mind.

**As built (Unit 5, `include/NiTCAD/physics/`):**
- `Semiconductor::create(SemiconductorParameters)` validates a parameter set once (`invalid_input` naming the parameter:
  non-finite values; non-positive permittivity, Eg0, Nc300, Nv300, mu_max, N_ref, alpha or lifetime; negative Varshni
  alpha/beta or mu_min; mu_min > mu_max), so the model functions rely on it. `silicon_parameters` holds the legacy `SILICON`
  values that the first diode uses; electron affinity, saturation velocity, Auger, band-gap narrowing, effective masses and the
  other legacy material sets come with the units that use them. `band_gap_eV`, `conduction_band_dos`, `valence_band_dos` and
  `intrinsic_density` are the legacy Varshni, (T/300)^1.5 and sqrt(Nc Nv) exp(−Eg/2kT) formulas, with kT from
  `base::thermal_voltage`. Successful validation does not allocate (messages are built only on failure).
- `check_temperature(m, T)` returns `invalid_input` unless T is finite and positive, Eg(T) > 0, and for both carriers
  mu_max(T) = mu_max (T/300)^T_exponent is finite and ≥ mu_min. Past either limit the formulas still evaluate but stop
  meaning anything: n_i exceeds sqrt(Nc Nv), or mobility rises with doping (silicon: holes from about 857 K, electrons from
  about 953 K; Eg reaches zero near 2998 K). The legacy checks neither. Unit 6 (or whichever layer takes T as input) calls it.
- `caughey_thomas_mobility(m, carrier, N_total, T)` and `scharfetter_lifetime(m, carrier, N_total)` take the total ionised
  impurity N_A + N_D. Both depend only on doping and temperature, so they return values without partials (item 2).
- `srh_recombination(n, p, EquilibriumProduct, n_ie, tau_n, tau_p)` returns R, dR/dn and dR/dp. It is `constexpr` and
  inline, as it runs per node and Newton iteration. The equilibrium product carries its own partials with respect to n and p,
  so it is the legacy `recombination_fd` form; with the Boltzmann product (n_ie², partials zero) it is
  `recombination_boltzmann` operation for operation. That is the R4 hook: Fermi–Dirac adds a producer of
  `EquilibriumProduct`, and no caller changes.
- `statistics.hpp` holds every place statistics enter: `boltzmann_density(n_ie, eta)` (density and d/deta from the reduced
  potential), `boltzmann_equilibrium_product(n_ie)`, and `boltzmann_neutral_equilibrium(C, n_ie)` (n, p and eta of a
  charge-neutral node, for ohmic contacts and the initial guess; majority carrier from the square root, minority from mass
  action, eta = asinh(C/2n_ie), as in the legacy `contact_value` and initial guess, but arranged so that no intermediate
  overflows or underflows where the result is representable, see OLD / NEW below). Assembly (Unit 7) must call these and not
  write a statistics formula itself. All of them are homogeneous in the concentrations, so the assembler may pass values
  already divided by Ns.
- Errors and preconditions: parameter sets and the temperature range are user input (`std::expected`, `create` and
  `check_temperature`); temperature, doping and n_ie in the once-per-problem functions are `NITCAD_EXPECTS` preconditions
  (the device or solve layer validates them first; the temperature precondition lives in one private helper,
  `src/physics/temperature.hpp`, beside the one mu_max(T) formula); the per-iteration kernels (`srh_recombination`,
  `boltzmann_density`) have no checks. Limit: `boltzmann_equilibrium_product` returns n_ie², which underflows for n_ie below
  about 1.5e-154 (cryogenic, or divided by Ns); Boltzmann results there are not meaningful anyway.
- Changed from the legacy (OLD / NEW / REASON):
  - OLD `N = max(N, 1)` in mobility and lifetime. NEW: N ≥ 0 is a precondition, with no clamp. REASON: the clamp changed
    mu(0) by about 1e-15 and silently accepted negative input (6.7).
  - OLD: `bernoulli` and `clip_700` were in `kernels.hpp` beside the recombination. NEW: not in physics; `boltzmann_density`
    does not clip. REASON: the Bernoulli function belongs to the Scharfetter–Gummel flux (Unit 7), and bounding a Newton
    overshoot belongs to the solver (6.7). exp overflows above eta ≈ 709.78.
  - OLD: kB in eV/K was the rounded 8.617333262e-5. NEW: `base::k_B_eV_per_K`, the exact ratio k_B/q. REASON: Unit 2
    already pins it. Effect on n_i(300 K): about 4e-10 relative.
  - OLD neutral equilibrium `0.5 (C + sqrt(C² + 4 n_ie²))`, `n_ie² / n`, `asinh(C / 2 n_ie)` (scaled units). NEW
    `0.5 |C| + hypot(0.5 C, n_ie)`, `n_ie (n_ie / n)`, and ln|C| − ln n_ie when C / 2n_ie overflows. REASON (PR #9 review):
    the old forms overflow for |C| near DBL_MAX and return p = 0 once n_ie² underflows. In the normal range the values
    agree to rounding.
- **For Unit 9:** the legacy J(0.5 V) fixture (`_build` in `tests/test_device1d_native_gates.py` [verified]) runs with
  `auger=True` and band-gap narrowing on. BGN is exactly zero there, because the doping (1e17) is below `bgn_N0` = 1.3e17.
  Auger is not in Unit 5. An estimate [derived, not run]: at 1e17, Auger adds about 0.3% to the bulk recombination rate, and
  bulk recombination is a small part of this short-base diode's current, so Auger's effect on J should be far below the 1%
  gate. Unit 9 must confirm this by measurement, or the owner must add Auger. (Measured at Unit 9: Auger moves J(0.5 V) by
  4.3e-7 relative, 6.2 "As built (Unit 9)".)

## 6. Cross-cutting design

### 6.1 de Mari scaling

The legacy solves in scaled variables. Its `Inputs` struct carries `VT, ni, Ns, LD, J0, R0, eps`, with a
doping scale `Ns` that defaults to `max(|doping|, ni)` unless overridden
(`inputs.hpp`; `Device1D` comment on `Ns` [verified]). The legacy `Device3D` derives its scaling from the
maximum |doping| of whatever array it was built with, so two devices covering different slices silently
disagreed on units until an `Ns_override` was added (`CLAUDE.md`, "Physics/model conventions" [verified]).

Proposal: scaling is an explicit `assemble` input derived once per problem, never recomputed per slice
or hidden inside a device.

**Legacy scaling definitions (V3, read from `core/src/device1d/inputs.cpp:41-65` and
`core/src/device1d/device1d.cpp:175-192`, 29 [verified]):**

| Quantity | Definition |
|---|---|
| V_T | k·T/q |
| ε | ε_r(node 0) · ε₀, in F/cm (other nodes enter through the relative permittivity `et`, below) |
| Ns (concentration scale) | the override if one is given (> 0), else max(max over nodes of \|doping\|, n_i), with n_i of node 0's material at T |
| L_D (length scale) | sqrt(ε · V_T / (q · Ns)) |
| D0_REF | 1.0 (`kD0Ref`, commented "kernels.D0_REF"); its unit is not stated in the file |
| J0 (current scale) | q · D0_REF · Ns / L_D |
| R0 (rate scale) | D0_REF · Ns / L_D² |

Scaled quantities (same files):

- Coordinates `xs = x / L_D`; edge length `h[k] = xs[k+1] − xs[k]`; node control volume
  `dV[i] = (h[i−1] + h[i]) / 2`, with `dV[0] = h[0]/2` and `dV[N−1] = h[N−2]/2`.
- Doping `C = doping / Ns`; per-node `nie_s = nie / Ns`, where `nie` includes band-gap narrowing when enabled.
- Unknowns per node, in this interleaved order: ψ (in units of V_T, so ψ = φ/V_T), n/Ns, p/Ns
  (`F[3i]`, `F[3i+1]`, `F[3i+2]` at `device1d.cpp:1129,1144,1164`).
- Equilibrium carriers `n = nie_s · exp(ψ + s)`, `p = nie_s · exp(−ψ − s)` (`s` is a band-offset gauge shift, zero for
  a single material). Ohmic contact: `n0 − p0 = C`, `n0·p0 = nie_s²`, `psi0 = V/V_T + ln(n0/nie_s) − s`
  (`device1d.cpp:322-340`).
- Poisson row: `et[e]·(ψ[i+1] − ψ[i])/h[e] − et[e−1]·(ψ[i] − ψ[i−1])/h[e−1] − dV[i]·(n − p − C) = 0`
  (`device1d.cpp:597-599`), with `et` the harmonic mean of the relative permittivity ratio ε_i/ε_0 on edges
  (`inputs.cpp:152-153`); `et = 1` for one material.
- Continuity rows: edge diffusivity `dn_edge = hmean(μ_n0) · V_T / D0_REF`, flux
  `Jn = (dn_edge/h) · (n[k+1]·B(δ) − n[k]·B(−δ))` with `B` the Bernoulli function, and the row
  `Jn[eR] − Jn[eL] − R_s · dV = 0`, with `R_s = R / R0` (`device1d.cpp:1007,1031,1095,1144`).
- Scaled currents are converted with `J0` (for example `device1d.cpp:1222-1223`).

**Findings relevant to later units:**

- The Newton update cap and tolerance are in scaled ψ: `max_dpsi = 5.0`, `tol_update = 1e-8`. In the
  equilibrium Boltzmann loop the correction is clipped to ±`max_dpsi` and the convergence test
  `max|Δψ| < tol_update` is then measured on the clipped value (`device1d.cpp:613-627`). The legacy lesson that
  convergence must be judged on the full correction, not the damped one, is recorded in the legacy
  `ARCHITECTURE.md` (M34) for the stiff paths [verified, quoted]; whether this equilibrium loop is exposed to it
  is not analysed here. Unit 8 must decide the criterion explicitly. (Decided at Unit 8: the full correction, see 6.2.)
- `tol_residual = 1e-10` is declared in the options (`device.py`); where the C++ solver uses it was not read.
- The legacy `Device3D` scales from whatever array it was built with (described earlier in this section); the explicit
  `Ns` input of the new design is the fix.
- The constants used by the C++ core (`core/include/tcad/physics/materials.hpp`: `kQ = 1.602176634e-19`,
  `kKB = 1.380649e-23`, `kEPS0 = 8.8541878128e-14`) are the same values as `pytcad/constants.py`.

**As built (Unit 7, `include/NiTCAD/assemble/`):**
- `make_scaling(device, Ns_override)` returns `Scaling` with the definitions above. The reference material is node 0's, as
  in the legacy. D0 is taken as 1 cm²/s: the legacy file states no unit, and cm²/s is the one that makes its scaled edge
  diffusivity μ V_T / D0 dimensionless. In D dimensions a control volume is divided by L_D^D and an edge's coupling area
  over length by L_D^(D−2), which reproduces the legacy 1D row exactly (tested against it written out by hand).
- `bernoulli` and `bernoulli_derivative` (inline). For |x| < 1e-2 they use a Taylor series; otherwise B = x/expm1(x) and
  B' = (1 − x − B)/expm1(x). The second form follows from B(−x) = B(x) + x, and it avoids dividing by x (MSVC Release
  reported C4723 for that). B is within 4 ulp and B' within 1e-13 of 50-digit references.
  OLD / NEW / REASON: the legacy clipped x to ±700, so B(−1000) was 700; NEW has no clip (expm1 gives B(x) = −x and
  B(x) = 0 at the extremes) because the clip was a silent clamp (6.7). The legacy switched to the series only below 1e-4,
  where its derivative formula lost about 4 digits; NEW switches at 1e-2 with a longer series.
- `sg_electron_flux` and `sg_hole_flux` return the legacy Jn = a(n₂B(δ) − n₁B(−δ)) and Jp = −a(p₂B(−δ) − p₁B(δ)), with
  δ = ψ₂ − ψ₁, together with exact partials with respect to both potentials and both densities. The edge factor a is the
  assembler's (continuity rows are Unit 9). Both fluxes vanish at Boltzmann equilibrium, and reversing an edge negates them
  exactly.
- `ohmic_contact_value(C, n_ie, V/V_T)` returns ψ0 = V/V_T + η and n0, p0 from `physics::boltzmann_neutral_equilibrium`.
- `EquilibriumPoisson::create(device, scaling)` precomputes the scaled volumes, couplings (times ε_r/ε_r,ref, the legacy
  `et`), doping, n_ie and contact values.
  - Contact rows are Dirichlet, ψ − ψ0.
  - The Jacobian pattern is the diagonal plus both directions of every edge, built once. Contact rows keep their
    off-diagonal entries as explicit zeros. `evaluate` only rewrites values, so a `LinearSolver` analyzes once over a
    Newton run (tested: three iterations, `analyses() == 1`).
  - `charge_neutral_potential()` is the legacy initial guess.
  - Heterojunctions (an edge between regions whose material parameters differ) are rejected as `invalid_input`: band
    offsets and permittivity steps are deferred. Disconnected regions may differ.
- Not in Unit 7: the Newton loop, its damping and `min_pivot_ratio` (Unit 8); continuity rows and bias (Unit 9). On an
  accepted device every connected part has a Dirichlet contact row (6.4, "As built"), so Unit 8 may set
  `min_pivot_ratio = 0` for this system.

### 6.2 Scaled Newton variables and tolerances

Newton works on scaled unknowns (ψ, then n and p in the form the legacy uses). Legacy defaults:
`tol_update = 1e-8` (maximum scaled update), `tol_residual = 1e-10`, `max_dpsi = 5.0` (damping cap on the
scaled potential update) (`pytcad/device.py`, solver options dataclass [verified]). Tolerances are stated in
scaled units. The public API stays in the cm/V convention; scaling is internal to `assemble`/`solve`.
Whether the new solver reproduces these exact defaults is to be confirmed unit by unit.

**As built (Unit 8, `include/NiTCAD/solve/`):**
- `newton_solve(system, x, NewtonOptions, LinearSolver&)` is a template over any system with `unknowns()`, `make_jacobian()`
  and `evaluate(x, F, J)`. Each iteration solves J dx = −F, clips every component of dx to ±`max_update`, and adds it to x.
  Defaults are the legacy ones: `max_iterations = 100`, `tol_update = 1e-8`, `max_update = 5`.
- **Convergence criterion (the 6.1 question, decided):** the largest component of the full correction, before clipping,
  must be below `tol_update`. The legacy tested the clipped correction. With a componentwise clip and `tol_update` <
  `max_update` the two tests agree, so the legacy equilibrium loop was not exposed, but the new test does not depend on that.
  A test pins the difference (a tolerance between the clipped and the full step does not stop the iteration).
- `tol_residual` is not used. The legacy 1D options declare 1e-10, but its 1D Newton never reads it; only `device2d.py` and
  `device3d.py` test a residual (`residual_ok`, 1e-7, in their own routine) [verified]. The final max|F| is reported.
- Errors: a non-finite residual or no convergence within `max_iterations` is `non_convergence` (context: iteration and
  last full correction); linear-solver errors pass through with the iteration prefixed to the message; bad options are
  `invalid_input`. x keeps the last iterate.
- `solve_equilibrium(device, EquilibriumOptions)` builds the scaling and `EquilibriumPoisson`, starts from the
  charge-neutral potential and returns `EquilibriumSolution` in V and cm⁻³ (6.2). `EquilibriumSolution` is plain data in
  `solve` until the results layer exists (Unit 10). The linear solver keeps its default configuration, including the
  pivot-ratio check: every equilibrium Poisson row has a Dirichlet anchor or a strong charge term, so there is no reason to
  turn it off (6.10; Unit 9 revisits it for continuity rows). Cancellation and progress are Unit 10.
- **Finding, built-in potential gate:** the legacy gate takes ψ(last) − ψ(first); both are Dirichlet contact values, so it
  checks the boundary values, not the solve. Unit 8 keeps it and adds the peak field of symmetric abrupt junctions against
  the depletion approximation with the 2 V_T correction (0.01–0.17% at 1e16–1e18, gate 0.5%; without the correction about 3%).
  Asymmetric junctions are outside that approximation: at 1e18/1e16, majority holes spill about 4 nm (a few Debye lengths)
  into the light side, about 1e11 cm⁻² against a depletion charge of 3e11 cm⁻², and the junction field is 39% above it.

**As built (Unit 9, `assemble/drift_diffusion.hpp`, `solve/bias.hpp`):**
- `assemble::DriftDiffusion` is the coupled system with unknowns (ψ, n, p) per node, interleaved 3i, 3i+1, 3i+2, as in the
  legacy. The rows are those of `Device1D::residual_jacobian` (baseline models), generalised to the mesh graph:
  - Poisson, as in `EquilibriumPoisson`, but with n and p as unknowns.
  - Electron and hole continuity: the sum of the SG fluxes over the node's edges, ∓ V·R.
  - SG edge factor: the harmonic mean of the two nodes' Caughey–Thomas mobilities, times V_T/D0, times the scaled
    coupling (the legacy `dn_edge/h`).
  - R: `physics::srh_recombination` on scaled densities, times Ns/R0 (SRH is homogeneous of degree one).
  - Contact rows: Dirichlet ψ − ψ0, n − n0, p − p0, at each contact's bias.
  - Jacobian pattern: a 3×3 block per node plus 5 entries per edge direction, built once.
- Jn and Jp are the conventional electron and hole currents along each edge. Their sum is divergence-free at every
  non-contact node, so a contact's terminal current is the total current on the edges leaving it, in A·cm^(D−3) (A/cm² in
  1D). The currents through all contacts sum to zero.
- Newton hooks (`update_size`, `apply_update`):
  - ψ is clipped to ±`max_update`; n and p are clamped to [0.1, 10]× their current value (legacy).
  - Convergence is judged on the full correction: max(|dψ|, |dn|/n, |dp|/p). The legacy measured the clamped one.
  - `NewtonReport` records the smallest pivot ratio.
- `DriftDiffusionModels{doping_mobility, srh}` are the legacy `Models` flags of the baseline diode. OLD / NEW / REASON:
  with `doping_mobility` off, the legacy used `mu_max` at 300 K whatever the temperature; NEW uses Caughey–Thomas at N = 0,
  so μ_max(T); REASON: the same model at all temperatures (identical at 300 K).
- `solve_bias(device, bias_V, options, initial)` starts from the given state (a sweep's previous point) or from the
  equilibrium solution, sets the contact nodes to their biased values (legacy) and returns `BiasSolution`: potential and
  densities, terminal currents, edge currents and the Newton report. These are plain data in `solve` until Unit 10.
- Not carried from the legacy solver: the line search and generation-strength stages (impact ionization and BTBT only),
  Robin contacts with surface recombination velocity, the energy-balance block, lagged field mobility, and the ln(n_ie)
  term in δ (band-gap narrowing is deferred, so n_ie is uniform within a material).
- **Auger (the Unit 5 question): not needed.** Measured with a temporary local patch adding the legacy Auger term (removed
  afterwards): J(0.5 V) on the gate fixture goes from 1.280076828e-2 to 1.280077381e-2 A/cm², 4.3e-7 relative.
- **Pivot ratio on reverse bias (the 6.10 gate): not flagged.** The smallest pivot ratio over −0.5, −2 and −8 V is 2.1e-3,
  against the 1e-11 threshold, so the default check stays on.
- **Finding, FD-Jacobian gate:** normalizing each column by its largest entry hides small terms. Dropping dR/dn from the
  Jacobian still passes the gate, because in a density column the SG coefficients (about D/h²) dwarf it. Unit 9 checks the
  recombination part separately (J(srh) − J(no srh) against FD of F(srh) − F(no srh), at a uniform state where every flux
  vanishes) and checks the 1D rows against the legacy formulas written out by hand on an asymmetric diode.
- **Finding, current continuity:** the total current is the same on every edge to 7.3e-8 (legacy gate 1e-6). The electron
  current alone, with SRH off, varies by up to 4.9e-6 along the device. That is rounding, not convergence: in the neutral
  n region the electron flux is the difference of one-sided terms 2.7e9 times larger.

### 6.3 Device description

A `device` is plain data: regions (geometry in the mesh's coordinates), doping per node or region,
material per region, and contacts (6.4). It contains no solver state and can be built without a solver
present. Per-region materials are in scope from the start of the data model, although heterojunction
physics is deferred.

### 6.4 Contact description

A contact is a named set of boundary nodes with a type. First unit: ohmic, with the charge-neutral
equilibrium boundary value. Schottky, gate and oxide-coupled contacts are deferred. Contact boundary
conditions are residual and Jacobian rows, so contact handling lives in `assemble`, and terminal current
is computed from fluxes in `assemble`/`solve`. It is not an `analysis` quantity.

**As built (Unit 6, `include/NiTCAD/device/`, sections 6.3 and 6.4):**
- `Device::create(DeviceDescription)` validates once and owns its parts: the mesh (moved in), a uniform lattice temperature
  (in the device, as in the legacy `Device1D(x, doping, T, …)`), named regions each holding a validated
  `physics::Semiconductor`, the region of every node, the donor and acceptor concentrations of every node, and the contacts.
  Accessors return spans; `net_doping`, `total_impurity` and `material` give per-node values.
- Doping is stored as N_D and N_A, not as net doping. Poisson needs N_D − N_A and the mobility and lifetime need N_D + N_A,
  so a compensated region is described exactly. The legacy took the net doping plus an optional `Ntotal` defaulting to
  |net| (`device.py`, `Device1D.__init__` [verified]); that default is the case where each node has one dopant type.
- `Contact` is a name, `ContactKind::ohmic` and strictly increasing node ids. Each node must be on a mesh boundary patch (6.4:
  boundary nodes) and in no other contact (two Dirichlet values for one node). The applied bias is not device data; it
  belongs to the solve layer.
- Every error is `invalid_input`, and the message and context index name the region, node or contact at fault: missing or
  repeated region names; the temperature, through `physics::check_temperature` for each region's material (the message
  names the region); per-node arrays of the wrong size, a region index out of range, a region with no nodes; doping not
  finite and ≥ 0; contact name, node and boundary errors.
- **Topology (the gate from 6.10):** every connected part of the mesh graph must contain a contact node, otherwise
  `invalid_input` ("floating region") with the lowest node of that part. Ohmic contact rows are Dirichlet rows for all
  three equations, so on an accepted device every connected part of every equation has one. Units 7 and 8 may then set
  `min_pivot_ratio = 0` for device solves, provided their contact rows really are Dirichlet. One contact is enough.
- Not carried from the legacy: `check_mesh` (it prints the worst spacing-to-Debye-length ratio and does not validate
  anything), the warning for doping above 1e19 cm⁻³ under Boltzmann statistics (there is no warning channel; Fermi–Dirac is
  deferred), and `graded_mesh` (not requested).

### 6.5 Input representation

The authoritative input is the in-memory C++ `device` description plus a solver-options struct. An
on-disk input format is **deferred** (R2) until concrete requirements exist. The input types must not depend on
any file format: no serialization types or format concepts inside them, and any future reader or writer depends on
them, never the reverse. The legacy wire format was a Python JSON `DeviceSpec` and is not carried over.

### 6.6 Result representation

`results` holds fields indexed by mesh node, terminal quantities per bias point, and a run record
(inputs hash/identity, options, per-iteration convergence history, converged flag). It is plain data.
An on-disk result format is **deferred** (R2), under the same independence rule as 6.5.

### 6.7 Error policy (R1, decided)

- **Recoverable errors** are returned as `std::expected<T, NiTCAD::base::Error>`: invalid input, a degenerate
  mesh, a singular or inaccurate linear solve, Newton non-convergence, cancellation. The caller decides what
  to do. Loud failure on invalid input is kept: nothing is silently clamped.
- **Exceptions** are used only for exceptional conditions and programmer errors (violated preconditions or
  invariants, allocation failure). Which mechanism enforces programmer-error preconditions (exception,
  assertion or terminate) is open (Q2).
- **No exceptions in hot numerical kernels or Newton iterations.** Model evaluation, assembly loops, the
  Newton loop and the linear-solver interface report failure by value. Where a dependency can throw, the
  layer that calls it converts the exception to an error value at its boundary.
- `Error` is a small typed value: a category (invalid input, degenerate mesh, singular system, inaccurate
  solve, non-convergence, cancelled), a message, and optional context. The category list is a proposal, to be
  confirmed at Unit 2.
- The legacy uses a typed C++ exception hierarchy (`core/include/tcad/base/errors.hpp` [verified]); the new
  design keeps its categories as values, not as exception types. Two legacy rules are kept: a clamp used
  during Newton overshoot must not be applied to the final converged value, and a linear solve whose accuracy
  check fails is an error. The legacy check was relative residual above 1e-6, for PARDISO only [verified, `CLAUDE.md`];
  the new one is the componentwise backward error, for every backend (6.10, "As built", and the review that motivated it).
- `std::expected` needs `<expected>` in the pinned MSVC toolset; verified in V1.
- **Precondition macro (Q2, decided):** one project macro, `NITCAD_EXPECTS(cond)`, always on in debug and release. On
  violation it logs file, line and condition, then fails fast (`std::abort` or `__fastfail`). It never throws.
  Preconditions are tested by testing the predicates, since Catch2 has no death tests. `__fastfail` vs `abort`
  behaviour under a debugger is not tested. Keep the macro thin: C++26 contracts may replace it.
- **Error shape (Q4, decided):** `enum class ErrorCode : std::uint8_t { invalid_input, degenerate_mesh,
  singular_system, inaccurate_solve, non_convergence, cancelled, resource_exhausted }` and
  `struct Error { ErrorCode code; std::string message; std::optional<ErrorContext> context; }`, where `ErrorContext`
  holds an optional index (node, edge, row or bias step) and an optional numeric value (residual or tolerance).
  There is no "internal error" code; internal bugs go through `NITCAD_EXPECTS`. The `std::string` allocates only on
  the error path, never on the success path.

### 6.8 Determinism and reproducibility

- Agreement between runs, machines and implementations is defined by **tolerances**, never by bit-identity
  digests. The legacy found that its bit-identity goldens pin one machine's floating-point summation order
  and must be regenerated per machine (`CLAUDE.md`, "Python/testing" [verified]); they are not carried over.
- Same input, same thread count, same build: results must be reproducible.
- A kernel may use threads only if its result is identical across thread counts (parallelize over the output
  index, never accumulate in thread order). Legacy default is 1 thread for that reason
  (`core/include/tcad/runtime/threads.hpp` [verified]).
- A multithreaded linear solver may differ at round-off between thread counts (legacy `CLAUDE.md` on MKL
  PARDISO [verified]); results then need tolerance-based comparison.

### 6.9 Threading, worker execution, cancellation, progress

- The numerical core is single-threaded by default (6.8).
- The application runs every solve on a **worker thread**, never on the UI thread. The worker owns the
  solver objects; the UI thread owns the windows.
- **Cancellation:** `solve` accepts a cancellation token (C++20 `std::stop_token` is the candidate) checked
  at safe points: between Newton iterations and between bias steps. A cancelled run returns a typed
  "cancelled" outcome and keeps the partial results. Cancellation is not possible inside a single linear
  factorization or solve. That latency is a stated limit.
- **Progress:** `solve` reports through a callback or observer interface (bias point, iteration,
  residual norm, convergence flag). The solver never calls UI code. The worker forwards progress to the UI
  thread by posting a Win32 message.
- This requires no extra layer: the token and observer types are defined in `solve`; the `app` layer owns the thread.

### 6.10 Linear solver interface (D3)

`linalg` exposes a small backend-neutral interface. Eigen SparseLU is the first backend. PARDISO and
iterative backends are deferred. The interface must be shaped so adding them needs no redesign of callers:

1. The input is a sparse matrix in a `linalg`-owned storage type plus right-hand sides. Eigen types are an
   implementation detail of the Eigen backend and do not appear in the interface.
2. Pattern analysis, numeric factorization and solve are separate steps, so a backend can reuse symbolic
   analysis across Newton iterations. The legacy `DirectSession` keeps its symbolic analysis across Newton
   iterations (`CLAUDE.md` [verified]). A backend that cannot reuse it simply repeats the analysis.
3. Backend choice and options (tolerance, iteration limit, preconditioner) are explicit
   configuration values, not environment variables. The legacy used environment variables such as
   `PYTCAD_LINSOLVE_BACKEND` (`CLAUDE.md` [verified]); they are not carried over.
4. Every solve returns `std::expected` with diagnostics (backward error, relative residual, refinement steps, and
   iteration count where it applies). Direct and iterative backends share one result shape. A solve that fails the
   accuracy check is an error even if the backend itself reported success (legacy: PARDISO perturbs tiny pivots
   instead of failing [verified]).
5. Thread count is part of backend configuration and defaults to 1 (6.8).

**Neutrality check (Q5, decided):** Unit 3 adds no second backend. Neutrality is proved mechanically and behaviourally:
a build test compiles the public `linalg` headers in a translation unit with no Eigen include path, so an Eigen type
leaking into the interface fails the build; and all Unit 3 tests are written against the interface only. Whether one
interface fits both direct and iterative backends stays unproven until one is requested.

**As built (Unit 3, `include/NiTCAD/linalg/`):**
- `SparseMatrix`: CSR with `std::int32_t` indices (Eigen's and LP64 PARDISO's index type), built from triplets.
  Duplicates are summed in input order; explicit zeros stay in the pattern; values can be rewritten in place, the pattern
  cannot. Bad dimensions or indices are `invalid_input`. It is `BasicSparseMatrix<Scalar>` for `Scalar` = `double`
  (`SparseMatrix`) or `std::complex<double>` (`ComplexSparseMatrix`), both compiled once in `sparse_matrix.cpp` (A5).
- `LinearSolver::create(SolverConfig)` validates `backend`, `threads` (Eigen SparseLU: exactly 1), `equilibrate`
  (default on), `max_backward_error` (default 1e-8), `max_refinement_steps` (default 1) and `min_pivot_ratio`
  (default 1e-11; must be 0 when not equilibrating).
- `analyze`, `factorize` and `solve` are separate. `factorize` re-analyzes whenever the pattern differs from the analyzed
  one, so a caller is always correct and the `analyses()` and `factorizations()` counters show the reuse. A failed analysis
  or factorization keeps nothing, so the next call analyzes afresh (legacy behaviour). With Eigen SparseLU the reusable part
  is only COLAMD and the elimination tree (2–6% of a 2D factorization, 0.1% of a 3D one, measured); the rest is redone
  every time because row pivoting depends on the values.
- Everything below lives in `LinearSolver`, so every backend gets it identically (linalg hardening, after the Unit 3 review):
  - **Equilibration:** rows, then columns, scaled by powers of two (exact). On a coupled system with rows scaled 1e±20
    and columns 1e±10 it took the forward error from 5.6e-4 to 4.3e-12.
  - **Singularity:** an empty row or column (analysis), an all-zero row or column, an exactly zero pivot, or a
    smallest-to-largest pivot ratio below `min_pivot_ratio` is `singular_system`, with the row or the original column as
    the index. The ratio catches a floating region, whose last pivot is at rounding level rather than zero. It is a
    heuristic. Measured on equilibrated systems: floating regions 1e-16 to 2.3e-14; well-anchored valid systems at
    least 2.7e-8, including graded meshes up to a 1e8 spacing ratio. Weakly anchored valid systems overlap the singular
    range: a region tied to the rest only by a zeroth-order term δ (relative to its couplings) gives about n·δ/2, and one
    joined by a single weak link w about w/2. With n = 1000, δ = 1e-14 or w = 1e-12 is rejected, although both solve
    with backward error 1e-16; an SRH-only floating body on a fine mesh can reach this range. A non-finite pivot
    (overflow during elimination) is `inaccurate_solve`.
  - **Acceptance:** the componentwise backward error max_i |b − Ax|_i / (|A||x| + |b|)_i (Oettli–Prager, with the
    Arioli–Demmel–Duff denominator in column-equilibrated units where that one is at rounding level) must not exceed
    `max_backward_error`, after up to `max_refinement_steps` of iterative refinement (stopping at the first step that
    does not reduce it, keeping the better iterate). It replaced `||Ax − b||₂/||b||₂ ≤ 1e-6`,
    which measured the condition number rather than the solve: it rejected a backward-stable 1D Laplacian with
    n = 600,000 and accepted a solve with 0.06% forward error on a badly scaled system. `||Ax − b||∞/||b||∞` is still
    reported.
  - Non-finite A or b is `invalid_input` and names the row. `std::bad_alloc` is caught at this boundary and becomes
    `resource_exhausted`. Shape mismatches, overlapping b and x, and use of a moved-from solver go through `NITCAD_EXPECTS`.
  - Measured cost of these checks: 1–5% of an Eigen factorization (2D 30k–120k unknowns, 3D 24k), within timing noise.
- Backend: Eigen SparseLU with COLAMD ordering, in `src/linalg/eigen_sparse_lu.*` (private), linked `PRIVATE`. It keeps a
  column-compressed copy of the pattern and a CSR-to-CSC position map, so each factorization only scatters values. It
  `static_assert`s Eigen 5.0.1, because its failure handling depends on that version's internals (message text, `info()`
  unset on allocation failure, the supernodal pivot storage, COLAMD checked only by `eigen_assert`). Eigen sizes COLAMD's
  workspace in 32-bit arithmetic and, when COLAMD fails, writes the permutation out of bounds before returning, so the
  backend refuses matrices whose workspace (about 2.2·nnz + 11·n words, computed in 64 bits) exceeds `INT32_MAX`
  (`resource_exhausted`; about 9.7e8 nonzeros) before calling it. It also checks the resulting ordering is a permutation,
  as a last line of defence.
- Not carried from the legacy `ReusableLU`/`DirectSession`: the subset-scatter and union-pattern re-analysis (only the
  deferred nonlocal models need them), the AMD ordering switch (measured pathological), and the size and instability
  fallbacks to SciPy (no SciPy here).

**Known limits, recorded for later units (from the Unit 3 review):**
- Eigen SparseLU does not scale to 3D: 81k unknowns (3D, 3 per node) took 70 s per factorization with 123× fill; 2D 270k
  took 6.8 s. A faster backend (PARDISO first; iterative or MUMPS for 3D, per the legacy measurements) is needed before 3D
  device work.
- Eigen stores L and U offsets in the 32-bit index type; extrapolated, that overflows near 5e5 3D unknowns, where a
  factorization would already take over 15 minutes. Fix with the next backend, or a 64-bit index in the private Eigen copy.
- About four copies of A exist during a factorization (caller, solver, CSC copy, Eigen's own), small next to L and U.
- The interface has one real right-hand side and an output-only `x`. Terminal admittances and sensitivities need several
  right-hand sides; AC small-signal solves a complex system (A5: the matrix type exists, a complex `LinearSolver` comes with
  the AC unit); iterative backends need an initial guess and the 3-unknowns-per-node block size.
- Pattern stability is the assembler's job: skipping a zero entry changes the pattern and forces a re-analysis. Unit 7
  should build the pattern once and assemble values in place; Newton tests should assert `analyses() == 1`.
- The pivot-ratio check cannot tell a floating region from a weakly anchored one (see Singularity above). Gates for later
  units: Unit 6 checks topology (every connected region of each equation has a Dirichlet row or a zeroth-order term),
  after which device solves may set `min_pivot_ratio = 0` (done: every connected part needs a contact, see 6.4 "As built"); Unit 9 must show that reverse-biased diode Jacobians (one-sided
  Scharfetter–Gummel links of order e^-40) are not flagged (done: smallest ratio 2.1e-3 at up to −8 V, 6.2 "As built (Unit 9)").
- Untested failure paths: `std::bad_alloc`, Eigen's out-of-memory messages, the COLAMD size guard (it needs about 1e9
  nonzeros), and a refinement step that increases the backward error.

### 6.11 Mesh generality (R3: unstructured meshes deferred, no later rewrite)

Unstructured meshes are not implemented initially. The mesh architecture must still not assume a tensor grid:

1. The public mesh interface is the graph: nodes (coordinates, control volume), edges (endpoint nodes, length,
   coupling measure) and boundary node sets. No layer above `mesh` uses (i, j, k) indexing, stride arithmetic
   or axis-aligned stencils.
2. Anything the legacy structured kernels compute per grid axis (for example `ii_grid.py` averages incident
   edges per axis [verified]) must be computed from edge direction data supplied by the mesh, not from grid axes.
3. Contacts and boundary conditions reference node (and, if needed, edge or face) sets by index, not index ranges.
4. The tensor grid is one producer of this graph. An unstructured producer later adds a constructor and its
   geometric-validity checks, nothing else. The legacy records one hazard to design against: its unstructured
   stencil gave negative dual-cell areas for clockwise triangles, a quirk pinned by a test and not fixed (legacy
   `ARCHITECTURE.md`, section 5.1, "OPEN DECISION from P2" [verified]).
5. How the coupling measure is defined on general meshes (for example a dual-cell construction) is not decided here.

**As built (Unit 4, `include/NiTCAD/mesh/`):**
- `Mesh` is the graph only: `points()` (x, y, z in cm, zero beyond the dimension), `volumes()`, `edges()` (`first < second`,
  `length`, `coupling_area`; a flux is area · (u_second − u_first) / length) and `boundary()` (named patches: strictly
  increasing node ids with a boundary face area per node). Edge directions come from the coordinates (item 2). Node ids are
  `std::int32_t`. In D dimensions volumes are in cm^D and areas in cm^(D−1) (1D per unit cross-section, 2D per unit depth).
- `Mesh::from_parts` is the single entry for every producer and validates the graph: dimension 1–3; finite coordinates;
  positive control volumes, edge lengths, coupling areas and boundary areas (`degenerate_mesh`, which also rejects the
  legacy negative-dual-area hazard of item 4); edge length equal to the endpoints' distance (relative 1e-12); endpoints in
  range, ordered, not repeated; patch names unique (`invalid_input`). The context index names the offending entry.
- `make_tensor_grid(x[, y[, z]])` is the first producer, following the legacy box method (`mesh2d.control_volume_widths`,
  `device2d.py:875–880`, `device3d.py:212–224`): widths are half of each adjacent spacing, a node's volume is the product of
  its widths, an edge's coupling area is the product of the other axes' widths (1 in 1D), and patches `x_min` … `z_max` carry
  the product of the face's widths. Node order (x fastest) is internal to the producer.
- Not carried from the legacy mesh modules: `graded_mesh` and `merge_mesh` (node placement, to come with the device unit if
  requested) and `check_mesh` (needs doping, so it belongs above `mesh`); the unstructured stencils (`stencil.cpp`) stay
  deferred (R3).

## 7. Language placement (D6: C and Fortran deferred)

| Component | Language | Reason |
|---|---|---|
| base, linalg, mesh, physics, device, results, assemble, solve, analysis, render, app | C++23 | primary language |
| Any later low-level kernel | C or Fortran | deferred; only for a measured, justified component named by the owner |

C and Fortran integration is deferred. If introduced later, it is used only for a measured and justified
component, behind a stable boundary (a C ABI), and the primary architecture stays C++23. Such a component must
be toolchain-compatible with D1 (MSVC); no Fortran compiler is chosen now. The legacy has no Fortran or C
source [verified: no `.f90/.c/.h` files in the legacy tree listing].

## 8. Repository layout

```text
ARCHITECTURE.md
DECISIONS.md
CMakeLists.txt              build; toolchain gate; vcpkg pin check           (Unit 1, exists)
CMakePresets.json           debug and release presets (Q3)                   (Unit 1, exists)
vcpkg.json                  manifest: baseline and version overrides (D8)    (Unit 1, exists)
.gitignore                  ignores build/                                   (Unit 1, exists)
tests/scaffold_test.cpp     Unit 1 scaffold test                             (Unit 1, exists)
include/NiTCAD/base/        constants.hpp, error.hpp, contract.hpp           (Unit 2, exists)
src/base/contract.cpp       NITCAD_EXPECTS failure path                      (Unit 2, exists)
tests/base/                 constants, error, contract tests and probe       (Unit 2, exists)
include/NiTCAD/linalg/      sparse_matrix.hpp, linear_solver.hpp             (Unit 3, exists)
src/linalg/                 sparse matrix, solver, private Eigen backend     (Unit 3, exists)
tests/linalg/               matrix, solver and header-boundary tests         (Unit 3, exists)
include/NiTCAD/mesh/        mesh.hpp, tensor_grid.hpp                        (Unit 4, exists)
src/mesh/                   graph validation, tensor-grid producer           (Unit 4, exists)
tests/mesh/                 graph validation and geometry gates              (Unit 4, exists)
include/NiTCAD/physics/     semiconductor, mobility, recombination, statistics (Unit 5, exists)
src/physics/                parameter validation, band quantities, models    (Unit 5, exists)
tests/physics/              published values, limits, FD derivative gates     (Unit 5, exists)
include/NiTCAD/device/      device.hpp, contact.hpp                          (Unit 6, exists)
src/device/                 description validation, topology check           (Unit 6, exists)
tests/device/               construction, validation, floating regions       (Unit 6, exists)
include/NiTCAD/assemble/    scaling, bernoulli, sg_flux, ohmic, equilibrium_poisson (Unit 7); drift_diffusion (Unit 9)
src/assemble/               scaling, scaled device (private), equilibrium Poisson, drift-diffusion (Units 7, 9)
tests/assemble/             Bernoulli/SG references, FD-Jacobian gate, reduction (Unit 7, exists)
include/NiTCAD/solve/       newton.hpp, equilibrium.hpp (Unit 8); bias.hpp (Unit 9)
src/solve/                  equilibrium and bias solves                      (Units 8, 9)
tests/solve/                Newton contract, equilibrium and bias diode gates, legacy graded_mesh port (Units 8, 9)
.github/workflows/ci.yml    CI: build and test Debug and Release per branch  (CI unit, exists)
include/NiTCAD/<layer>/...  public headers per layer                         (created per unit)
src/<layer>/...             implementations                                  (created per unit)
tests/<layer>/...           Catch2 tests mirroring src/                      (created per unit)
app/                        Win32 entry point                                (created at L8 only)
```

The layer names are decided (D5) and the namespace and include root are decided (Q1). Not part of the structure until
the owner approves them in a unit request: `cmake/` (toolchain modules), `benchmarks/` (measured performance cases). A unit adds its headers, sources and tests together.
`tests/scaffold_test.cpp` sits outside `tests/<layer>/` because Unit 1 has no layer.

## 9. Build, test and warnings

- **Build system:** CMake. The C++23 requirement is enforced by the build, not trusted to defaults. The language
  mode flag is passed explicitly (see Toolchain below), and configuration fails if a probe using the C++23 features
  the code actually uses (for example `<expected>`) does not compile in that mode.
- **Toolchain (D1, Q6, Q9: decided):** MSVC. The language mode is **`/std:c++23preview`, set explicitly**.
  `/std:c++latest` is not used, and `CMAKE_CXX_STANDARD 23` is not relied on, because CMake 4.3.1 maps it to
  `-std:c++latest` (rolling, `_MSVC_LANG` 202400). The build must:
  1. pass `/std:c++23preview` explicitly;
  2. fail configuration unless `_MSVC_LANG == 202302L` and the feature macros the code needs are defined
     (`__cpp_lib_expected`, `__cpp_lib_jthread`, and others as units adopt them). Verified in a scratch
     CMake project (V4);
  3. fail configuration if the compiler silently ignores or rejects the C++23 flag. **Verified here:** an
     unknown `/std:` value produces only warning D9002 and the build continues; `/WX` and `/we9002` do **not**
     turn it into an error. `/options:strict` does (error D8043). So the requirement is met by `/options:strict` plus
     the `_MSVC_LANG` check, not by `/WX` alone.
- **Toolset pin (Q9, decided):** the verified reference toolset is **MSVC 14.51** (toolset directory
  `14.51.36231`, compiler 19.51.36260), with Windows SDK `10.0.26100.0`, CMake `4.3.1-msvc1` and Ninja `1.13.2`
  as found on this machine. CI pins exact MSVC toolset, Windows SDK, CMake and Ninja versions, uses Visual Studio
  **Build Tools** rather than the full IDE, and must reproduce these versions (not yet verified, Q10). Local
  Visual Studio Community may be used for individual development subject to Microsoft's licence terms.
  The owner's earlier Q9 target of **14.50 LTS** was **not adopted**: 14.50 is not installed here and was not
  tested, so it is not recorded as supported. 14.51 is a regular release with a shorter support window (section 12), so
  moving to an LTS toolset is a planned later step (Q11).
- **Dependencies (D8, Q7, Q8: decided):** CMake with vcpkg in manifest mode, a pinned baseline, explicit version
  constraints and overrides where required.
  - Baseline: the full commit `fbb0f7bb200b07a9eb9081c7a3cf51d1aa1c51a1`, the one tested in V2. vcpkg HEAD is not
    tracked. The baseline changes only through an explicit dependency-update change followed by full validation.
  - Eigen **5.0.1** and Catch2 **3.16.0** come from that baseline. Both built and passed a probe test with
    `/std:c++23preview` (V4).
  - Triplet: **`x64-windows-static-md`** for development, CI and release, unless a concrete requirement justifies
    another.
  - Verified locally (V2): a manifest with `builtin-baseline`, `version>=` and `overrides` resolved with the Visual
    Studio-bundled vcpkg. The CI cache strategy is not specified (Q10).
- **Flags in effect (Unit 1, `CMakeLists.txt`):** `/std:c++23preview /options:strict /EHsc /utf-8 /permissive-
  /fp:precise /W4 /WX /external:W0`, definitions `UNICODE _UNICODE NOMINMAX WIN32_LEAN_AND_MEAN`, and the dynamic
  CRT (`MultiThreaded` or `MultiThreadedDebug` + `DLL`). CMake's default `/W3` is removed so it does not conflict. No
  `/arch:` flag, so no AVX2 assumption. N1–N4 (14.2) are provisional defaults.
- **Warning policy:** first-party code is compiled at `/W4` with warnings treated as errors (`/WX`); imported
  headers are system headers whose warnings are silenced (`/external:W0`). `/options:strict` is part of the policy.
- **Vcpkg tool pin (Q3):** `CMakeLists.txt` requires `VCPKG_ROOT` to be a git checkout whose `HEAD` equals the
  manifest's `builtin-baseline`, so the tool and the port versions cannot drift apart. Eigen 5.0.1 and Catch2 3.16.0 are
  pinned by `overrides` in `vcpkg.json`.
- **CI (Q10, `.github/workflows/ci.yml`):** runs on every branch push, every pull request and on demand; one
  job per preset (Debug and Release) on `windows-2025`. It checks out the vcpkg tool at the manifest baseline,
  installs CMake 4.3.1 and Ninja 1.13.2 from downloads checked against pinned SHA-256 values, selects the toolset with
  `vcvarsall.bat x64 -vcvars_ver=14.51`, and **fails** unless `cl` is 19.51.36260, the Windows SDK is 10.0.26100.0,
  CMake is 4.3.1 and Ninja is 1.13.2. Actions are pinned by commit (`actions/checkout` v5.1.0, `actions/cache` v5.1.0).
  vcpkg binaries use the `files` source in a directory saved with `actions/cache`, keyed on triplet, toolset and
  `vcpkg.json`. The CMake checksum is the one Kitware publishes (verified); Ninja publishes none, so its value was
  computed locally from the v1.13.2 asset. **First run (windows-2025, Visual Studio 2026 Enterprise 18.10.2) verified:** `-vcvars_ver=14.51` selects
  compiler 19.51.36260, the SDK is 10.0.26100.0, CMake is 4.3.1 and Ninja is 1.13.2. It also showed that the Visual Studio
  developer environment overwrites `VCPKG_ROOT` with its bundled vcpkg; the workflow now restores the pinned checkout after
  `vcvarsall.bat`, and so must anyone building from a Visual Studio developer prompt. **Cache verified on a re-run of the same
  workflow:** the first run missed and saved; the second restored all 4 packages in 850 ms (vcpkg install step 228 ms). Whole-job time
  barely changed (about 2m20s to 2m35s, Debug and Release), so the cache is not the dominant cost.
- **Test levels** per component, from the legacy practice: analytic and limiting cases → published-value
  regression → Jacobian vs finite differences → dimensional-reduction identity → convergence and mesh
  independence → benchmark. A compile-only check is never the sole gate for a numerical unit.

## 10. Validation criteria pinned from the legacy

These are the acceptance criteria the new units should meet. Source and verification status are stated for each.

| Criterion | Value | Source / status |
|---|---|---|
| Constants | Values as in `pytcad/constants.py`: q = 1.602176634e-19 C, k = 1.380649e-23 J/K, ε₀ = 8.8541878128e-14 F/cm, ħ = 1.054571817e-34 J·s, m₀ = 9.1093837015e-31 kg | [verified] file contents. The file does not name the edition; the numbers are the CODATA 2018 values [verified against NIST, V7]. Unit 2 pins CODATA 2018 in `constants.hpp`. |
| Thermal voltage | V_T(300 K) = 0.025852 V | [derived] 0.0258519998 V from the legacy k and q |
| n_i(300 K), Si | ≈ 1.0674e10 cm⁻³ | [derived] hand-computed from legacy `materials.py` (Eg0 = 1.17 eV, α = 4.73e-4 eV/K, β = 636 K, Nc300 = 2.86e19, Nv300 = 3.10e19, k = 8.617333262e-5 eV/K): Eg(300) = 1.12452 eV, n_i = 1.06738e10. **Not a legacy test assertion.** The legacy test only requires 9e9 < n_i < 1.6e10 (`tests/test_model_benchmarks.py`, `test_ni_300k_within_accepted_band` [verified]). Proposed new gate: 1.0674e10 with relative tolerance 1e-4, which holds if the same parameters are used. |
| Built-in potential | within 2e-3 V of V_T ln(Nd·Na/n_i²) | [verified] `tests/test_device1d_native_gates.py` (1e17/1e17 fixture). It measures ψ(last) − ψ(first), two Dirichlet contact values, so it does not test the solve; Unit 8 adds the peak-field gate (6.2, "As built"). |
| Diode J(0.5 V) | 1.280e-2 A/cm² within 1% | [verified] `tests/test_device1d_native_gates.py::test_g2_forward_current_matches_documented_value`, fixture `_diode()` in that file. The same test records the short-base analytic value 1.321e-2. Fixture parameters are to be extracted when that unit is requested. |
| Diode law | J/J_ideal in (0.85, 1.15) at 0.5 V; ideality within ±0.02 of 1 for V ≥ 0.3 V | [verified] `tests/test_validation.py::test_ideal_diode_law`. That test uses Caughey–Thomas mobility and has SRH on by default (`Models` dataclass, `device.py`: `doping_mobility=True`, `srh=True` [verified]); this is why both come before the bias solve. |
| Jacobian vs finite differences | worst per-column-normalized error ≤ 5e-5 over ≥ 80 columns | [verified] `tests/test_m13_solver.py` ("house FD-Jacobian gate (<= 5e-5, >= 80 columns)"). Normalization [verified at Unit 7, `_jacobian_probe`]: at a state perturbed from the solution (ψ by 0.02 normal noise, n and p by 1%), central differences with step 1e-7·max(\|u\|, 1) for randomly chosen columns; per column, max over rows of \|FD − J\| divided by the column's largest \|J\| + 1e-30. |
| Dimensional reduction | a transversely uniform 2D solution is y-independent to atol 1e-9 | [verified] `tests/test_validation_2d.py`. Legacy also reports 1.11e-16 V for 3D → 2D (`examples/05_3d_reduces_to_2d.py`, quoted in legacy `ARCHITECTURE.md` [unverified: not re-run]). New tolerances are set per unit. |
| Mesh geometry | total volume/length matches the domain to 1e-14 (2D) / 1e-10 (3D) | [verified] `tests/test_validation_2d.py`, `tests/test_validation_3d.py` |

## 11. Proposed build order (each unit = one branch)

Unit 1 was named by the owner on 2026-10-02 and has no open prerequisites (V1–V4 are done; V3 is consumed at Unit 7). The
order is by dependency, so each unit is
testable on arrival. Units 4–9 are dimension-generic from the start (D4) and use no tensor-grid assumptions (6.11).

| # | Unit | Layer | Legacy reference | Gate |
|---|---|---|---|---|
| 1 | Build scaffold: CMake, C++23 gate, vcpkg manifest with pinned baseline, Catch2 v3 harness | — | `core/CMakeLists.txt` (reference for options only) | configure fails on a toolset without the needed C++23 features; dependencies resolve from the pinned manifest; one trivial test runs |
| 2 | Constants, units, `Error` type, `NITCAD_EXPECTS` (**done, on `main`**) | base | `constants.py`, `core/include/tcad/base/errors.hpp` | CODATA 2018 values; V_T(300 K) = 0.025852 V |
| 3 | Sparse matrix, backend-neutral solver interface (6.10), Eigen SparseLU backend (**done, on `main`, with the hardening**) | linalg | `linsolve.py`, `core/src/solver/direct_lu.cpp` | known systems (e.g. analytic tridiagonal); accuracy check (backward error after hardening); singular system returns an error; symbolic-reuse path exercised; header-boundary check, no second backend (Q5) |
| 4 | Generic node/edge/control-volume mesh; tensor-grid constructors for D = 1, 2, 3 (**done, on `main`**) | mesh | `mesh.py`, `mesh2d.py`, `mesh3d.py`, `core/include/tcad/mesh/stencil.hpp` | total volume matches domain (section 10); positive dual volumes; consistent edge geometry across D = 1, 2, 3 on a uniform grid; public interface contains no (i, j, k) indexing |
| 5 | Si material parameters, Caughey–Thomas mobility, SRH recombination, Boltzmann statistics (value + partials) (**done, on `main`**) | physics | `materials.py`, `core/include/tcad/physics/` | n_i(300 K) ≈ 1.0674e10; published mobility values (legacy `test_caughey_thomas_matches_published_silicon_values`); derivatives vs finite differences; SRH vanishes at equilibrium |
| 6 | Device description and ohmic contact data (**done, on `main`**) | device | `core/include/tcad/device1d/inputs.hpp` (reference only) | construction and validation; invalid input returns an error |
| 7 | Scaling (V3), SG/Bernoulli flux, equilibrium Poisson residual + Jacobian, ohmic boundary (**done, on `main`**) | assemble | `core/src/device1d/device1d.cpp`, `core/src/device1d/inputs.cpp` | Bernoulli limits and symmetry; FD-Jacobian gate (section 10) |
| 8 | Newton solver (scaled variables) and equilibrium solve (**done, on `main`**) | solve | `device.py` options, `device1d.cpp` | built-in potential within 2e-3 V; bulk neutrality; convergence; non-convergence returns an error value |
| 9 | Electron/hole continuity assembly and bias solve; first end-to-end gate (**done on branch `solve/bias-continuity`**) | assemble, solve | `device1d.cpp`, `tests/test_validation.py`, `tests/test_device1d_native_gates.py` | J(0.5 V) = 1.280e-2 A/cm² ± 1%; ideal-diode law; current continuity; mesh independence; uniform 2D/3D reproduces 1D to a tolerance set at this unit |
| 10 | Result representation, cancellation and progress | results, solve | — (new) | cancellation returns partial results; progress is monotonic |
| 11+ | Everything else: Fermi–Dirac, Auger, band-gap narrowing, field mobility, heterojunctions, impact ionization, BTBT, transient, AC, thermal, process, unstructured meshes, PARDISO/iterative backends, file formats, analysis, render, app | deferred | per the audit | per unit, when the owner requests |

Units 1–9 are the smallest end-to-end vertical slice: a validated drift-diffusion diode. Everything past
Unit 10 is deferred until sequenced. Units 10 onward are not started.

## 12. Legacy facts: verified, derived and unverified

**Verified against the legacy files**
- Dimension-generic 2D/3D kernels exist: `ii_grid.py`, `dg_grid.py`, `thermal_grid.py`, `btbt_grid.py`,
  `hydro_grid.py`. Each states it is one kernel for Device2D and Device3D, and they are 2D/3D, not 1D.
  `dg_grid.py` states that only its Λ rows are shared; the classical Poisson row stays written per device.
  So the earlier wording "legacy 1D/2D/3D are three codebases" was too strong: the core Poisson/DD assembly
  is per-device while several optional physics kernels are generic. The unstructured DD modules also exist
  (`unstructured_dd.py`, `unstructured_dd3d.py`; file listing).
- "98% in sparse LU": origin is legacy `ARCHITECTURE.md` line 268 ("24^3 equilibrium solve spends 98% wall
  time in `_superlu.gssv`; assembly 0.011 s of 0.494 s"). That is SciPy SuperLU on a **3D equilibrium**
  solve of 0.494 s total. It is not Eigen and not drift-diffusion. It supports "the direct-solver
  algorithm dominated wall time there", not a claim about any new solver.
- PARDISO vs Eigen: `core/include/tcad/solver/pardiso_lu.hpp` lines 11–13: "Measured 2026-09-29 on B3's 2D
  MOSFET Jacobian (72,912 unknowns): refactor+solve 218 ms (1 thread) / 99 ms (10 threads) vs Eigen SparseLU
  1,456 ms." This is a recorded comment, not re-run here. `benchmarks/BASELINE.md` lists B3 at 11,640 DOF for
  its quick size, so the 72,912-unknown figure is a different size.
- Legacy C++ standard in its core CMake is C++20 (`core/CMakeLists.txt`, line 39). MKL is loaded at
  runtime from `mkl_rt` (`CLAUDE.md`).
- Defaults `doping_mobility=True`, `srh=True`; Newton defaults `tol_update=1e-8`, `tol_residual=1e-10`,
  `max_dpsi=5.0` (`pytcad/device.py`).
- `RecombinationResult{R, dRdn, dRdp}` in `core/include/tcad/physics/kernels.hpp` (value plus exact partials).
- Every test criterion in section 10 marked [verified].

- **V3 scaling:** the definitions and scaled-variable forms in 6.1 were read from `inputs.cpp` and `device1d.cpp`.
- **V1 (MSVC, run here 2026-10-02):** toolset and `/std` behaviour as in section 9; with `/std:c++23preview` and
  `/std:c++latest` the probe defines `__cpp_lib_expected` 202211, `__cpp_lib_mdspan` 202207, `__cpp_lib_print` 202406,
  `__cpp_lib_jthread` 201911, `__cpp_multidimensional_subscript` 202211, `__cpp_explicit_this_parameter` 202110,
  and compiles and runs a program using `std::expected`, `std::mdspan`, explicit object parameters, `std::println`
  and `m[i, j]`. With `/std:c++20` the same program fails (no `<expected>` contents, no explicit object parameters).
  `<stop_token>`/`std::jthread` is available in C++20 mode.
- **V2 (vcpkg, run here):** Visual Studio-bundled vcpkg (`2026-07-27-98d7cb0…`) with `VCPKG_ROOT` pointed at it;
  a manifest pinned to baseline `fbb0f7bb200b07a9eb9081c7a3cf51d1aa1c51a1` (the head of microsoft/vcpkg on
  2026-10-02, used for the test, not chosen) resolved `catch2` 3.16.0 and `eigen3` 5.0.1. An override to `eigen3`
  3.4.0 with a `version>=` constraint on Catch2 also resolved (Eigen 3.4.0, Catch2 3.16.0). Both builds passed the
  probe test under `-std:c++latest`. Catch2 was consumed as the `Catch2::Catch2WithMain` target.

- **V4 (run 2026-10-02):** a scratch CMake project with the explicit `/std:c++23preview` flag, a configure-time
  `try_compile` gate checking `_MSVC_LANG == 202302L`, `__cpp_lib_expected` and `__cpp_lib_jthread`, Eigen 5.0.1 and
  Catch2 3.16.0 from the pinned baseline: configured, built with `/W4`, and the Eigen SparseLU and Catch2 test passed.
  `/WX` and `/we9002` leave an ignored `/std:c++23` as a built executable; `/options:strict` yields error D8043.
- **MSVC support lifecycle** (Microsoft C++ blog, "New Release Cadence and Support Lifecycle for MSVC Build Tools",
  fetched 2026-10-02 [verified: the page's stated policy]): a new MSVC release every six months (May and
  November) with 9 months of support; each second November release is LTS with 3 years of support. 14.50
  (November 2025) is the LTS release; 14.51 (May 2026) is a regular release; 14.52 (November 2026) is the next LTS.
  The page gives no end dates. The end dates used in this document (about February 2027 for 14.51, about
  November 2028 for 14.50) are **[derived]** from those periods.

**Derived here, not a legacy assertion**
- V_T(300 K) = 0.0258519998 V; n_i(300 K) = 1.06738e10 cm⁻³ (inputs in section 10).

**Unverified: do not rely on these**
- The unit of `D0_REF`: unstated in the legacy; Unit 7 takes it as 1 cm²/s (6.1, "As built"). Read since: the 5e-5
  Jacobian gate's normalization (Unit 7, section 10) and where `tol_residual` is used (Unit 8, 6.2).
- clang-cl and MinGW claims about C++23 (not tested; only MSVC 14.51 was probed). **Anything about MSVC 14.50**
  (not installed here, not tested), including its C++23 feature support. A search result states that 14.52 will
  include a `/std:c++23` switch and be C++23 complete; this was not verified and 14.52 is not released.
  A CI runner's Visual Studio Build Tools version, and whether CI can reproduce the pinned toolset.
- Fortran toolchain claims: that MSVC has no Fortran compiler; which Fortran compilers (Intel ifx, LLVM Flang,
  gfortran) can interoperate with which C++ ABI on Windows; CMake generator support for each.
- The MKL licence terms. Whether a CI runner can reproduce the local vcpkg build (cache behaviour).
- Whether Direct3D 12 / Direct2D headers work under a MinGW toolchain.

## 13. Risks

1. **C++23 on Windows.** Verified: this toolset has no stable `/std:c++23`, and an unknown `/std` value silently
   falls back to C++14 with only a warning that `/WX` does not promote. Mitigation (decided, Q6): explicit
   `/std:c++23preview`, `/options:strict`, and a configure-time check of `_MSVC_LANG` and the feature macros.
   Residual risk: the preview mode may change between toolsets, and 14.51 has a short support window (Q11).
2. **Fortran coupling to the C++ toolchain.** If Fortran is ever used, its compiler and runtime must be
   ABI-compatible with the C++ compiler, which constrains D1. Mitigation: D6 defers Fortran, and any later
   Fortran sits behind a C ABI in a separate library.
3. **Linear solver choice.** The legacy measured a large PARDISO speedup over Eigen SparseLU at one size
   (section 12), but PARDISO is proprietary and depends on MKL. Mitigation: baseline Eigen behind a small
   `linalg` interface; revisit with a new-repo benchmark.
4. **Generic mesh design cost.** Larger initial design effort than a per-dimension port. The legacy has the
   reduction tests and generic kernels as evidence it is workable (section 12).
5. **Unvalidated legacy physics.** Hydrodynamic, MC implant, TED and self-heating are simplifications with
   gaps disclosed in the audit. They stay deferred.
6. **Direct3D 12 verbosity.** Large and low-level, but isolated at L7/L8 and built last.
7. **Working tree under OneDrive.** The repository is under `OneDrive\Desktop`, so build output in `build/` is inside a
   synced folder. Not measured; Unit 1's first configure (vcpkg install) took about two minutes. If sync causes locks
   or slowness, `binaryDir` in `CMakePresets.json` can move out of the repository.

## 14. Decisions

Recorded by the owner. D7 of an earlier draft (the branch name for this file) was not a decision and is
withdrawn; numbering is otherwise kept stable. This proposal lives on `architecture/architecture-proposal`.

### 14.1 Decided

| # | Decision |
|---|---|
| D2 | Catch2 v3 is the C++ test framework. Its availability through vcpkg is checked under V2; the decision stands. |
| D4 | Dimension-generic architecture for 1D, 2D and 3D. |
| D5 | Root namespace `NiTCAD`. Layers: `base`, `linalg`, `mesh`, `physics`, `device`, `assemble`, `solve`, `results`, `analysis`, `render`, `app`. Nesting and include root: Q1. |
| R1 | `std::expected` for recoverable errors; exceptions only for exceptional or programmer-error cases; no exceptions in hot numerical kernels or Newton iterations (6.7). |
| Q1 | Layers are nested namespaces `NiTCAD::<layer>`; the include root is `include/NiTCAD/<layer>/`, spelled exactly like the namespace. Settle casing before any file exists, because Git `core.ignorecase` can hide case-only renames. |
| Q2 | `NITCAD_EXPECTS(cond)` fail-fast precondition macro (6.7). |
| Q3 | A standalone vcpkg tool checkout pinned to a recorded commit (not the Visual Studio-bundled copy); `CMakePresets.json` is part of Unit 1; the CI workflow is not part of Unit 1 and gets its own owner-named unit. |
| Q4 | `ErrorCode` enum (seven values) plus `Error` struct (6.7). |
| Q5 | No second backend in Unit 3; neutrality proved by a header-boundary build check and interface-only tests (6.10). |
| Q10 | CI design direction: a GitHub-hosted Windows runner, toolset selected with `vcvarsall.bat x64 -vcvars_ver=14.51`, exact versions asserted (compiler 19.51.36260, SDK 10.0.26100.0, CMake 4.3.1, Ninja 1.13.2), vcpkg `files` binary cache saved with `actions/cache` and keyed on baseline, triplet, toolset and manifest hash. **Design only, not verified:** runner contents, `vcvars_ver` selection on a runner and the cache behaviour were not run. Implemented in `.github/workflows/ci.yml` on branch `architecture/ci-workflow` (the CI unit, not Unit 1); **not yet run on a runner**, first run done, see section 9. |
| Q11 | Move to an LTS toolset when MSVC 14.52 LTS is released (expected November 2026), no later than 60 days after, and before 14.51 support ends (about February 2027, derived). Re-run V1 and V4, run the full suite, update the CI assertions and cache key, on one dedicated branch. |
| R4 | Boltzmann carrier statistics initially; Fermi–Dirac and other statistics deferred until requested (section 5, item 7). |
| D1 | MSVC is the Windows C++23 toolchain (C++ compiler only; Fortran is D6). Toolset pin per Q9. |
| Q6 | `/std:c++23preview` set explicitly; `/std:c++latest` not used. The build verifies `_MSVC_LANG` and the C++23 feature macros and fails configuration if the C++23 flag is ignored or invalid (9). |
| Q7 | The exact full vcpkg baseline commit tested in V2 (`fbb0f7bb200b07a9eb9081c7a3cf51d1aa1c51a1`); vcpkg HEAD is not tracked. Eigen 5.0.1. The baseline changes only through an explicit dependency-update change followed by full validation. |
| Q8 | `x64-windows-static-md` is the initial and standard triplet for development, CI and release unless a concrete requirement justifies another. |
| Q9 | Verified reference toolset MSVC 14.51 (`14.51.36231`), with exact MSVC toolset, Windows SDK, CMake and Ninja versions pinned for CI, Visual Studio Build Tools on CI, and local Community use subject to Microsoft's licence terms. The 14.50 LTS target was not adopted (untested); an LTS move is Q11. |
| D8 | CMake + vcpkg manifest mode with a pinned baseline, explicit version constraints/overrides, and a reproducible CI dependency/cache strategy. Baseline, Eigen version and triplet are decided (Q7, Q8); the CI cache strategy is open (Q10). |
| A5 | **AC representation and index width** (decided 2026-10-02, owner: "Do A5 now"). AC small-signal solves the complex system (J + iωC) x = b, as the legacy code does (`ac.py:187`, `ac2d.py:317`, `ac3d.py:256`: `spsolve(J0 + 1j*omega*Cmat, b)`); the real 2N block form is not used (4× the nonzeros and a departure from legacy). `SparseMatrix` is therefore `BasicSparseMatrix<Scalar>` with `Scalar` ∈ {`double`, `std::complex<double>`}, done while only linalg and its tests depend on it; `SparseMatrix` stays the name of the real matrix, so no caller changes. Assembly (Unit 7) should build J and, later, C on one shared pattern so the AC matrix reuses it. Complex solving (a complex `LinearSolver`, its acceptance check and backends) is part of the AC unit. Index width stays `std::int32_t`: it is the index of Eigen SparseLU and LP64 PARDISO and allows 2³¹ − 1 nonzeros in A, beyond the reach of a direct factorization; the backend-internal factor limits are in 6.10 "Known limits". |

### 14.2 Provisionally decided (stand unless verification or Unit work disproves them)

| # | Decision | Provisional because |
|---|---|---|
| D3 | Eigen SparseLU is the initial sparse direct solver behind a backend-neutral interface that allows PARDISO and iterative backends later without redesign (6.10). | The interface shape was built in Unit 3 (6.10, "As built") and awaits the owner's review of that unit. Fit for an iterative backend is unproven until one is requested. |
| N1 | Dynamic CRT: `CMAKE_MSVC_RUNTIME_LIBRARY` = `MultiThreaded$<$<CONFIG:Debug>:Debug>DLL`, matching `x64-windows-static-md`. | Advisor suggestion; the owner has not ruled on N1–N4 explicitly. |
| N2 | `/fp:precise` and no `/arch:` flag (no AVX2 assumption). | A neutral default chosen by me for Unit 1; the advisor gave no recommendation. Affects tolerance-based reproducibility (6.8). |
| N3 | `/EHsc` (Eigen and the standard library can throw). | Advisor suggestion. |
| N4 | `UNICODE`, `_UNICODE`, `NOMINMAX`, `WIN32_LEAN_AND_MEAN` defined and `/utf-8` set. | Advisor suggestion. |
| N5 | CI builds and tests both Debug and Release on every branch push and every pull request. | My default for the CI unit; the owner asked for CI on every feature branch but has not ruled on the configurations. |

### 14.3 Deferred

| # | Item | Condition to revisit |
|---|---|---|
| D6 | C and Fortran integration | A measured, justified component is named by the owner. The primary architecture stays C++23. |
| R2 | External input/output file formats | Concrete requirements exist. Internal input/result types stay independent of any format. |
| R3 | Unstructured-mesh implementation | The owner requests it. The mesh design must already allow it without a rewrite (6.11). |
| — | Fermi–Dirac and other statistics | Requested (R4). |
| — | PARDISO and iterative solver backends | Requested; the interface already accommodates them (6.10). |

### 14.4 Verification results (run 2026-10-02; probes were built outside the repository)

| # | Result | Status |
|---|---|---|
| V1 | Toolset: VS 18 Community, MSVC dir `14.51.36231`, compiler 19.51.36260. Required features are present (`std::expected`, `std::mdspan`, explicit object parameters, multidimensional subscript, `std::println`, `<stop_token>`/`jthread`) under `/std:c++23preview` and `/std:c++latest`. **No stable `/std:c++23`**; it is silently ignored (D9002, falls back to C++14). CMake 4.3.1 emits `-std:c++latest` for `CXX_STANDARD 23`. Details in section 12. | Verified. Q6 decided. |
| V2 | The Visual Studio-bundled vcpkg works in manifest mode with a pinned baseline, `version>=` constraints and overrides. Eigen (5.0.1 at the test baseline, 3.4.0 by override) and Catch2 3.16.0 built and passed a probe test under `-std:c++latest` (superseded by V4, which uses the decided `/std:c++23preview`). Not verified: CI reproduction and cache behaviour, a standalone vcpkg checkout, whether the redirect of vcpkg's build/download/package directories out of Program Files was necessary (done as a precaution). | Verified locally. Q7, Q8 decided. CI reproduction open (Q10). |
| V4 | Re-check under the decided mode: Eigen 5.0.1 and Catch2 3.16.0 from the pinned baseline build and pass under explicit `/std:c++23preview` with a configure-time `_MSVC_LANG` gate (passed). `/WX` and `/we9002` do not stop an ignored `/std:c++23`; `/options:strict` does (D8043). MSVC 14.50 was **not** tested (not installed). | Verified on 14.51 only. |
| V5 | Q3 prerequisite and Unit 1 run: a standalone vcpkg checkout at the baseline commit (partial clone, `--filter=blob:none`, in `C:\Users\disha\vcpkg`, bootstrapped; tool release 2026-09-26) resolved the manifest with the Eigen 5.0.1 and Catch2 3.16.0 overrides. The Unit 1 configure, build and test passed in Debug and Release. Negative checks: `/std:c++20` fails the gate (`#error`, C1189); an ignored `/std:c++23` fails it (D8043); a `VCPKG_ROOT` that is not a git checkout at the pinned commit fails configuration. **Not tested:** an unset `VCPKG_ROOT`, CI, any other machine. | Verified locally. |
| V6 | CI workflow: the YAML parses and defines the intended triggers, matrix and steps; the pinned CMake download matches Kitware's published SHA-256 and Ninja runs as 1.13.2 locally. First CI run (push and pull request, run 36939739285): all version assertions passed; the job failed on the `VCPKG_ROOT` override described in section 9 (fixed in `a2b4b6a`). After the fix both Debug and Release passed on push and pull request, and the vcpkg binary cache restored on a re-run. | Verified on the hosted runner. |
| V7 | Unit 2 (`core/base`): the legacy constants are the CODATA 2018 values: eps0 8.8541878128(13)e-12 F/m, m_e 9.1093837015(28)e-31 kg and the exact hbar, per NIST (`physics.nist.gov/cuu/pdf/wall_2018.pdf`; CODATA 2022 changed eps0 and m_e at about 1e-9, so pinning 2018 is what preserves legacy results). `FAST_FAIL_FATAL_APP_EXIT` = 7 in Windows SDK 10.0.26100.0 `winnt.h`. Debug and Release: 12 Catch2 test cases pass, no warnings. A wrong eps0 (CODATA 2022) fails the constants test. | Verified locally. |
| V8 | Unit 3 (`core/linalg`): Debug and Release build with no warnings and pass 23 Catch2 test cases in `nitcad_linalg_test`. Neutrality check bites: a public header including `<Eigen/Core>` (C1083), including `<eigen3/Eigen/Core>` (C1189 on `EIGEN_WORLD_VERSION`) and Eigen linked `PUBLIC` (C1189 on `__has_include`) each fail the build; checked by temporary edits, then reverted. Eigen 5.0.1 SparseLU does not set `info()` when it cannot allocate its working memory but always sets `lastErrorMessage()`; the backend tests the message first. **Not tested:** `std::bad_alloc` and Eigen's out-of-memory path (`resource_exhausted`), more than 2³¹ − 1 nonzeros. | Verified locally. |
| V9 | linalg hardening (`core/linalg-hardening`): Debug and Release build with no warnings; `nitcad_linalg_test` has 35 test cases, all pass (0.45 s in Release). Regression tests reproduce each review finding: the n = 600,000 Laplacian is accepted (backward error ≤ 1e-15, ‖r‖/‖b‖ > 1e-6); the badly scaled coupled system is rejected without equilibration (backward error > 1e-6) and accepted with it (forward error ≤ 1e-10); floating regions are `singular_system` (pivot ratio < 1e-13), and the column reported for one floating region among contacted ones lies in it; a moved-from solver is neither analyzed nor factorized. Found while building: the Arioli–Demmel–Duff fallback with plain ‖x‖∞ fired on every row of a column-scaled system and hid a bad solve; it now uses column-equilibrated units. Pivot-ratio calibration and the cost of the checks are in 6.10. **Not tested:** `std::bad_alloc`, Eigen's out-of-memory path, a COLAMD failure. | Verified locally. |
| V10 | Follow-up review of the hardening (same branch): the pivot ratio read from Eigen equals the hand value (2/3 for diag(1, 3)); the singular index lies in the offending 4-node block at columns 3, 417 and 995, for both exact and rounding-level singularity; weakly anchored valid regions overlap the singular range (recorded in 6.10, pinned by a test); Eigen's COLAMD workspace is sized in 32 bits and its failure path writes out of bounds (`Eigen_Colamd.h:267`, `Ordering.h:140`), now guarded before the call; NaN pivots were ignored by the ratio scan, now `inaccurate_solve` (tested with an overflowing 2×2); refinement keeps the better iterate. Debug and Release: 39 test cases pass (0.45 s Release), no warnings; includes two solvers on two threads matching a serial run bit for bit. | Verified locally. |
| V11 | Unit 4 (`core/mesh`): Debug and Release build with no warnings; `nitcad_mesh_test` has 16 test cases, all pass. Gates of section 11 on graded axes (1e3 spacing ratio): total volume equals the domain to 1e-14 relative in 1D and 2D and 1e-13 in 3D (legacy: 1e-14 2D, 1e-10 3D absolute); all volumes, lengths and areas positive; 2D and 3D grids reduce to the 1D one when grouped by x (node volumes and x-edge areas equal the 1D values times the transverse area, to 1e-14); the box method reproduces the discrete Gauss identity for linear fields in 1D, 2D and 3D to 1e-12, and that check detects a 1% error in a single coupling area; each boundary face's areas sum to the face's measure. The public interface has no grid indices: every test reads only the graph. | Verified locally. |
| V12 | Unit 5 (`physics/silicon-models`): Debug and Release build with no warnings; `nitcad_physics_test` has 33 test cases, all pass. Gates of section 11: n_i(300 K) = 1.06738e10 (1e-4 gate; 1e-13 against a 40-digit reference computed independently from the legacy formulas); Caughey–Thomas passes the legacy published-value test (mu_n(0) = 1360, mu_n(1e18) = 263 within 25% of 300) and matches 40-digit values at 1e15–1e20 cm⁻³ and at 400 K to 1e-13; SRH is zero at equilibrium (exactly for an exactly representable np = n_ie², else at the rounding level of np) and gives dp/tau_p and dn/tau_n in low injection; dR/dn and dR/dp agree with central differences to 6.4e-11 worst (gate 1e-8), including the chain through a carrier-dependent equilibrium product. After the PR #9 review: `check_temperature` accepts silicon at 1–850 K and rejects T ≤ 0 or non-finite, Eg(T) ≤ 0 (including exactly 0) and mu_max(T) below mu_min (holes at 860 K, electrons at 1000 K) or infinite; neutral equilibrium stays finite and exact at |C| = DBL_MAX, with n_ie = 1e-170 (p = 1e-300 where the legacy gives 0) and where C / 2n_ie overflows; successful validation makes no heap allocation (counted through a replaced global operator new). Mutation checks, each failing exactly one test case: dropping the dE/dn term; alpha_n 0.91 → 0.90 (the legacy 25% mobility gate alone would not catch it; the pinned values do); each of the three legacy neutral-equilibrium forms; the old kT spelling (tested at temperatures where it rounds differently, 22–95.5 K); building the message on the success path; removing the mobility rule; accepting Eg = 0. | Verified locally. |
| V13 | Unit 6 (`device/device-description`): Debug and Release build with no warnings; `nitcad_device_test` has 14 test cases, all pass. Valid 1D and 2D silicon diodes (contacts taken from the `x_min`/`x_max` patches), a compensated node and a two-region device are read back unchanged. Each validation rule rejects its case with `invalid_input` and the documented index. Topology: a device with no contact, either half of two interleaved disconnected chains, and an isolated node are rejected, with the lowest node of the floating part as the index; one contact per part, or one contact spanning both parts, is accepted. Mutation checks: removing the topology check, the boundary rule, the shared-node rule, the unused-region rule, or the temperature check each fails the suite. | Verified locally. |
| V14 | Unit 7 (`assemble/equilibrium-poisson`): Debug and Release build with no warnings; `nitcad_assemble_test` has 16 test cases, all pass. Scaling of a 1e17 silicon diode matches 50-digit values (L_D = 1.29288e-6 cm, J0 = 12392.28 A/cm², R0 = 5.98249e28 cm⁻³s⁻¹) to 1e-14. Bernoulli: within 4 ulp (B) and 1e-13 (B') of 50-digit references from −1000 to 709; B(−x) = B(x) + x and B'(x) + B'(−x) = −1 hold to rounding; continuous across the series switch; no clipping. SG fluxes vanish at equilibrium to 1e-13 of their one-sided terms, reduce to diffusion and to upwinded drift, change sign exactly under reversal, and their partials match finite differences to 1e-8. Equilibrium Poisson: section 10 FD-Jacobian gate on every column (91 in 1D, 364 in 2D, 348 in 3D) at a perturbed state, worst 2.7e-9 (gate 5e-5); the 1D residual equals the legacy row written out by hand to 1e-13; a y-uniform 2D residual equals the 1D one times the scaled transverse width to 1e-12; contact rows are exact Dirichlet rows; three Newton iterations reuse one analysis. Mutation checks: dropping the charge derivative, a coupling scale of L_D^(D−1), edge terms in contact rows, a flipped hole-flux sign, and no heterojunction check each fail the suite. | Verified locally. |
| V15 | Unit 8 (`solve/equilibrium-newton`): Debug and Release build with no warnings; `nitcad_solve_test` has 12 test cases, all pass. Newton on one-unknown systems: quadratic convergence, one analysis per run, clipping with the last iterate kept, the full-correction criterion, option validation, non-finite residual and singular Jacobian errors. Equilibrium of a silicon diode on the legacy fixture mesh (spacing 1e-8 cm at the junction growing by 1.2 to 1e-6 cm, 241 nodes): converged in 7 iterations, the last corrections shrinking quadratically, final max|F| ≤ 1e-10, one analysis; V_bi within 2e-3 V at 300 K and 400 K; |p − n + C| ≤ 1e-6 |C| more than 0.4 µm from the junction; total charge ≤ 1e-8 of the depleted charge; n p = n_i² to 1e-13; peak field of symmetric junctions within 0.5% of the corrected depletion approximation (measured 0.01–0.17%); a y-uniform 2D and 3D device reproduce 1D to 1.1e-16 V (legacy gate 1e-9); an undoped device gives zero potential in one iteration; too few iterations, a bad Ns and a bad linear configuration return errors. Mutation checks: the clipped-correction criterion, no clipping, no finite-residual check and unscaled output potential each fail the suite. | Verified locally. |
| V16 | Unit 9 (`solve/bias-continuity`): Debug and Release build with no warnings; `nitcad_assemble_test` has 22 test cases and `nitcad_solve_test` 23, all pass. Fixtures on the legacy meshes, from a test port of `graded_mesh` that matches the legacy output (250, 262 and 450 nodes; nodes to 1e-12). **J(0.5 V) = 1.28008e-2 A/cm²** (gate 1.280e-2 ± 1%; 0.006% off) in 10 iterations, without Auger; Kirchhoff to 8.8e-9 and the total current equal on every edge to 7.3e-8 (legacy spread 1e-6); ideal-diode law: J/J_ideal(0.5 V) = 0.969 (gate 0.85–1.15) and ideality within 0.02 of 1 for every step from 0.3 to 0.6 V; least-squares ideality over 0.3–0.7 V 1.004 (gate ±0.05); a 2× finer mesh moves J(0.5 V) by 0.005% (gate 3%); y-uniform 2D and 3D devices reproduce the 1D current to 1.3e-15; reverse bias: J < 0 at −0.5, −2, −8 V with exponent 0.85 (gate 0.5–1.5), smallest pivot ratio 2.1e-3 (not flagged); zero bias is the equilibrium solution in one iteration; one analysis per Newton run, quadratic finish. Coupled FD-Jacobian gate (legacy probe at +0.3 V, every column; 183 in 1D, plus 2D and 3D): worst 3.3e-9; recombination part on its own: 7.6e-7 (gate 1e-5); the 1D rows equal the legacy formulas written out on a 1e18/1e16 diode to 1e-12; a y-uniform 2D residual equals the 1D one times the scaled width to 1e-12. Auger measured with a temporary patch: 4.3e-7 of J(0.5 V). Mutation checks: dropping dR/dn from the Jacobian, flipping the recombination sign in the hole row, flipping the terminal-current sign, an arithmetic instead of harmonic mobility mean, and no density clamp each fail the suite (the first and fourth only after the two targeted tests were added). | Verified locally. |
| V3 | Scaling definitions, scaled variables and Newton tolerances read from `inputs.cpp` and `device1d.cpp` and recorded in 6.1. One open question remains for Unit 8 (convergence criterion on the clipped correction). | Verified |

Verifications due at their own unit: none left. The 5e-5 Jacobian gate's normalization was read at Unit 7 (section 10), and
the `_diode()` fixtures at Unit 9 (`graded_mesh` with h_min 1e-8, h_max 1e-6, ratio 1.15 or 1.12; 1e17/1e17; CT, SRH, and
Auger and BGN that do not matter, 6.2).

### 14.5 Remaining architectural questions

| # | Question | Needed by |
|---|---|---|
| N6 | Licence of NiTCAD. Eigen is MPL-2.0, which affects static linking of release binaries. The repository `LICENSE` file was not checked | before release |
| N7 | Whether clang-format, `/analyze` or clang-tidy are wanted at all (feature freeze) | owner |

Nothing in this table blocks Units 1–3.

**Approval workflow:** the owner names each unit. Each unit is built on its own branch, tested, committed and reported, and is
not merged without the owner's instruction.
