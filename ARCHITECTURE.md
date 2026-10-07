# NiTCAD — Architecture

**CURRENT ARCHITECTURE, updated for Units 1–12.** Units 1–10 established the validated drift-diffusion foundation;
Unit 11 added band-gap narrowing and Auger, and Unit 12 added gate contacts on a lumped oxide (the MOS capacitor and the
MOSFET). Both are merged on `main` (Unit 12 as `1de5967`). The architecture below keeps the implemented layers and
numerical contracts intact while adding a long-term north-star for a general multiphysics semiconductor equation engine.
Sections distinguish implemented design from target architecture: the north-star material (2.1, 2.2, 6.12, 11.1, 11.2,
14.6) is future target architecture, not a new layer and not a claim about the current implementation. The dependency
architecture of sections 3 and 4 is authoritative.

The owner accepted the original architectural decisions D1–D6, D8 and R1–R5. Section 14 records their status.
This document must continue to be rewritten to describe what actually exists as units land; the north-star sections are
design direction, not permission to build features without an owner-named unit.

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
| Platform | Windows only; x64 target |
| Languages | C++23 (primary); C and Fortran only where technically justified |
| C++ standard | C++23, enforced by the build system, not by compiler defaults |
| Desktop layer | Native Windows desktop via Win32 API + Direct3D 12 + Direct2D; x64 target. No Qt, VTK, Python GUI or other GUI framework |
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
7. **NiTCAD is a general equation engine, not a collection of device-specific solvers.** A diode, MOS capacitor,
   MOSFET, heterojunction device, power device or future multiphysics device must be assembled from the same
   concepts: regions, interfaces, contacts, fields, equations, constitutive models and solver infrastructure.
8. **Physics is composable.** Adding a model should add or activate equations/terms and their derivatives; it should
   not require a new device solver or duplicate Newton, linear-algebra, mesh or result code.
9. **Interfaces are first-class physical objects.** Region-to-region and region-to-contact interfaces carry their
   own boundary conditions, displacement conditions, offsets, charge, recombination, transport or other interface
   models as required by the requested physics.
10. **Fields and equations are explicit concepts even when their current implementation is specialized.** The present
    code may expose `n`, `p` and `psi` through specialized result/system types, but future units should move toward
    reusable field and equation descriptors without breaking the existing numerical contracts.
11. **Numerical robustness has priority over feature count.** Every new physical capability needs convergence,
    conservation, Jacobian and dimensional/mesh gates appropriate to the equations it introduces.
12. **The north-star is additive.** The target architecture may introduce equation graphs, adaptive meshes, automatic
    differentiation and multiple nonlinear/linear solver strategies, but none is treated as implemented until an
    owner-named unit builds and verifies it.

### 2.1 North-star architecture: the equation graph

The long-term abstraction is an **equation graph** spanning `device`, `physics` and `assemble`. It is a conceptual
architecture, not a new dependency layer. The graph connects:

```text
Device description
  Regions ───────────────┐
  Interfaces ────────────┤
  Contacts / boundaries ─┤
                         ↓
                     Fields
          (unknown, fixed, derived)
                         ↓
                    Equations
      ┌──────────────────┼──────────────────┐
      ↓                  ↓                  ↓
   local terms       edge terms       interface/boundary terms
      └──────────────────┼──────────────────┘
                         ↓
                 Global F(x) and J(x)
                         ↓
                  Nonlinear solver
                         ↓
                   Linear solver
```

The key rule is that a new device type should normally be a new **composition of existing concepts**, not a new
solver implementation.

### 2.2 First-class concepts

| Concept | Long-term role | Current status |
|---|---|---|
| Region | geometry, material, doping and applicable equations | Implemented in `device`; a region is of a semiconductor or, since Unit 15b, of an insulator (Poisson only) |
| Interface | shared boundary between regions with interface laws | Semiconductor heterointerfaces since Unit 15 (band offsets through the band shift, permittivity steps, thermionic emission); semiconductor-insulator interfaces since Unit 15b (fixed charge, traps, surface recombination); the Unit 12 gate is a lumped oxide on a boundary patch |
| Contact / boundary | named boundary set plus physical boundary condition | Ohmic contacts, gate contacts on a lumped oxide (Unit 12) and electrodes on a meshed insulator (Unit 15b) implemented; Schottky contacts deferred |
| Field | named physical quantity with location, units and applicability | Present as specialized solution/result data; generic field abstraction is target |
| Equation | residual/Jacobian contribution associated with fields and domains | Present through specialized assembler systems; generic equation objects are target |
| Model | pure constitutive law with exact derivatives | Implemented for current local physics |
| Mesh | geometric graph supporting 1D/2D/3D and future unstructured producers | Implemented graph abstraction; unstructured/adaptive mesh deferred |
| Solver | nonlinear/linear algorithms independent of device type | Newton + backend-neutral linear solve implemented; more strategies deferred |
| Result | fields, terminals, convergence, provenance and diagnostics | Implemented plain-data results; richer scientific data model is target |

This table is a design direction. It does not add code to the current dependency graph until a unit requests it.

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

The **equation graph is not an additional layer**. It is a cross-layer design concept: `device` describes where
equations apply, `physics` provides constitutive models and derivatives, and `assemble` turns the active equations
and boundary/interface laws into the global residual and Jacobian. This preserves the existing one-way dependency graph.

## 4. Layer responsibilities

| Layer | Owns | Does NOT own |
|---|---|---|
| base | physical constants, unit convention, error types | physics |
| linalg | sparse matrix storage, backend-neutral solver interface (6.10), Eigen SparseLU backend (D3), residual check | meshes, Newton |
| mesh | node/edge/control-volume graph, per-edge length and coupling area, per-node volume; tensor-grid constructors for 1D/2D/3D; no tensor-grid indexing in its public interface (6.11) | doping, materials, linear algebra |
| physics | material parameters and constitutive models; node/edge/interface-local laws and exact derivatives | meshes, matrices, bias, solving |
| device | regions, material assignment, doping, interfaces and contact/boundary descriptions (data only) | solving, assembling |
| results | field/terminal data, convergence history, diagnostics and run record | solving |
| assemble | equation contributions, scaling, Scharfetter–Gummel flux (and its extensions), contact/interface residuals, F(x) and J(x), gathering model inputs | Newton iteration, UI |
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
5. **Interface-local and boundary-local models** may consume values gathered from the adjacent regions or boundary
   state, but they still return only their physical value and exact partial derivatives. The model does not own mesh
   traversal or Newton iteration.
6. **Nonlocal models** (path-integrated band-to-band tunnelling, nonlocal ionization) need paths across the
   mesh, so they live in `assemble`, not `physics`. They are deferred.
7. A new physics feature must be expressible as one or more reusable models/equation contributions. It must not
   create a device-specific nonlinear solver, duplicate contact handling, or duplicate linear-algebra infrastructure.
8. The first physics unit contained only what the first diode needed; later units may extend the model library without
   changing existing callers.
9. **Carrier statistics (R4):** Boltzmann statistics initially; Fermi–Dirac for parabolic bands since Unit 14 (an
   addition beside Boltzmann, below); others are deferred. Statistics stay isolated inside `physics`, so adding
   Fermi–Dirac was an addition, not a signature break for callers.
   Evidence from the legacy: `recombination_fd()` additionally takes the equilibrium product `np_eq` and its
   partials `dnpq_dn`, `dnpq_dp`, which `recombination_boltzmann()` does not
   (`core/include/tcad/physics/kernels.hpp` [verified]). Unit 5 was designed with that in mind.

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

**As built (Unit 11, band-gap narrowing and Auger; owner request after the Unit 10 review):**
- `SemiconductorParameters` gains `auger{Cn, Cp}` and `bandgap_narrowing{E0_eV, N0}`, with the legacy silicon values
  (2.8e-31 and 9.9e-32 cm⁶/s; 6.92 meV and 1.3e17 cm⁻³). Validation: the Auger coefficients and E0 must be ≥ 0, and N0 > 0.
  Setting them to zero switches the model off through the material.
- `bandgap_narrowing.hpp`: `bandgap_narrowing_eV(m, N_total)`, the legacy Slotboom form E0[ln(N/N0) + sqrt(ln²(N/N0) + ½) −
  sqrt(½)] above N0 and 0 at or below it (so there is no step at N0); `effective_intrinsic_density(m, N_total, T)` =
  n_i exp(ΔEg/2kT). OLD / NEW / REASON: the legacy clamped N to at least 1 cm⁻³; NEW requires N ≥ 0; REASON: as for
  mobility (Unit 5).
- `auger_recombination(n, p, EquilibriumProduct, Cn, Cp)` returns (Cn n + Cp p)(np − E) with exact partials through E (the
  legacy form). It is cubic, so it is called with physical densities, not scaled ones.
- `assemble::PhysicsModels{doping_mobility, srh, auger, bgn}` replaces `DriftDiffusionModels`; all four are on by default,
  as in the legacy. Both assemblers take it:
  - With `bgn`, every node's n_ie is the effective one. That changes the statistics, the SRH and Auger equilibrium
    product, the ohmic contact values and the charge-neutral guess.
  - The SG driving term becomes δn = Δψ + Δln n_ie for electrons and δp = Δψ − Δln n_ie for holes (the legacy `delta`,
    `delta_p`), so equilibrium carries no current where n_ie varies.
  - `EquilibriumOptions` gains `models` (only `bgn` matters there), and the run record lists `models.auger` and
    `models.bgn`.
- Effect on earlier gates: none at 1e17 (BGN is zero below 1.3e17, and Auger moves J(0.5 V) by 4.3e-7). The Unit 8
  peak-field gate at 1e18 now uses V_bi from the effective n_ie of each side, and the Unit 9 "no recombination" check now
  switches off both SRH and Auger.
- **Finding, BGN regimes:** BGN multiplies the diffusion current by exp(ΔEg/kT), measured 8.5908 against 8.5932 for
  1e19/1e19 at 0.5 V with SRH off. With SRH the low-bias current of such a junction is depletion-region recombination,
  which scales with n_ie, not n_ie²: measured 3.14 at 0.3 V, against √8.59 = 2.93.
- **Finding, current resolution in heavily doped devices:** extracted currents below about 1e-8 A/cm² are rounding noise
  in a 1e19 device. The majority one-sided fluxes reach 5e6 A/cm² even on 1e-6 cm cells, and the total current varies by
  5e-8 to 1e-7 A/cm² from edge to edge. This is the same cancellation as the Unit 9 electron-current floor and is
  inherent to the (ψ, n, p) unknowns the legacy also uses. Low-current analyses of heavily doped devices (leakage,
  sub-threshold) will need either a quasi-Fermi-potential formulation or current extraction from the minority side; this
  is not done here.

**As built (Unit 14, Fermi–Dirac statistics; owner request, roadmap item 14):**
- `physics/fermi_dirac.hpp`: `fermi_half(η)` returns F₁/₂ and F₋₁/₂ = dF₁/₂/dη (normalised to e^η as η → −∞);
  `log_degeneracy(η)` returns ln γ = ln F₁/₂ − η and its slope; `inverse_fermi_half(ν)`. Everything is evaluated
  through ln γ, so nothing overflows or underflows where the result is representable:
  - η ≤ −2: the exact series, summed for ln γ and its slope directly (no cancellation; at most 22 terms);
  - −2 < η ≤ 40: quintic Hermite interpolation of ln γ on a table of step 0.025;
  - η > 40: the Sommerfeld expansion with ten terms (residual 3.7e-18 at 40, falling like e^−η).
  The table holds ln γ and its first two derivatives at 1681 nodes. They come from 32-point Gauss–Legendre quadratures of
  F₁/₂, F₋₁/₂ and F₋₃/₂: t = s² on [0, 1], then panels of width ≤ 2 up to max(η, 0) + 40, summed with Neumaier
  compensation (naive sums left 4e-15). It is built once per process, on first use, as a function-local static (C++
  makes that initialisation thread-safe). It is immutable, a cached constant, and is the one exception to item 1's "no
  global state". Measured against 40-digit mpmath references (polylogarithm for η < 0, edge-subdivided quadrature
  above, checked against each other and against F₁/₂(0) = (1 − 2^−½) ζ(3/2)):
  - F₁/₂ is within 1.5e-15 everywhere.
  - F₋₁/₂ is within 1.3e-13 up to η = 15 and 3e-12 up to 40 (the slope F₋₁/₂/F₁/₂ falls to 0.04 there), and 1e-14
    above.
  F₋₁/₂ is the derivative of the function evaluated (the interpolant's own), so Jacobians built from it match finite
  differences of the residual. The inverse is safeguarded Newton on ln F₁/₂ = ln ν (private `increasing_root.hpp`,
  shared with the neutral equilibrium), 2–5 evaluations, and recovers η to 1e-14 relative from −700 to 1000.
- `statistics.hpp` adds the Fermi–Dirac producers beside the Boltzmann ones (the R4 hook), in the same gauge: η is
  measured from the intrinsic level, and with g = ln(N/n_ie) for the carrier's band,
  n = N F₁/₂(η − g) = n_ie e^η γ(η − g), which tends to n_ie e^η. With band-gap narrowing n_ie is the effective one,
  and g_n + g_p is the narrowed gap over kT. The producers:
  - `fermi_dirac_density(n_ie, g, η)`: the density and dn/dη, from the correctly rounded e^η times γ.
  - `Degeneracy` and `fermi_dirac_degeneracy(n_ie, g, density)`: ln γ at the density's own reduced energy, and
    d ln γ / d density. This is the legacy L and w of `fd_node_factors`.
  - `fermi_dirac_equilibrium_product(n_ie, Degeneracy_n, Degeneracy_p)`: n_ie² γ_n γ_p with its partials (legacy
    `npq_args`).
  - `fermi_dirac_neutral_equilibrium(C, n_ie, g_n, g_p)`: safeguarded Newton on ln(majority) = ln(|C| + minority),
    from the Boltzmann root; the majority is then |C| + minority, so neutrality holds to rounding.
  The degeneracy is a per-iteration kernel with no checks: a non-finite density gives NaN, so a diverging Newton
  iterate surfaces as a non-finite residual (6.7). Against 40-digit roots for silicon at 300 K, the neutral equilibrium
  is within 1e-13 in η and 1e-12 in n and p at ±1e20, 1e19, −3e19, 1e17 and 0. At 1e20 the Fermi level is 2.43 kT
  inside the band (legacy G7(a): > 2).
- OLD / NEW / REASON:
  - OLD (`fermi.py`, `fermi.hpp`): cubic Hermite tables of ln F₁/₂ and ln F₋₁/₂ on [−10, 40], step 0.005, with
    ln F₋₁/₂'s slope by an 8th-order finite difference, so F₋₁/₂ was not exactly the derivative of F₁/₂. Both
    functions refused outside [−40, 40], and the equilibrium solve clamped η at 40.
  - NEW: one quintic table of ln γ, with the series below it and Sommerfeld above it. Every η is valid; nothing is
    clamped or refused.
  - REASON: the derivative Newton uses is then exactly that of the residual; a parabolic band is an approximation
    near 5e21 cm⁻³, not a wall; and the bound on Newton overshoot belongs to the solver (6.7).

  - OLD (`fd_node_factors`): for η ≤ −30, L = w = 0 exactly, to keep `fd=True` bit-identical to Boltzmann there;
    densities were clamped at 1e-300.
  - NEW: no threshold. Below density/N = 1e-6, ln γ = −log1p(a v + b v²), with v = density/N, a = 2^−3/2 and
    b = ¼ − 3^−3/2. This avoids an inversion, its error is below 1e-18, and it is analytic through v = 0, so it also
    accepts a density at or just below zero (as a finite-difference probe of a minority density produces).
  - REASON: the threshold made w jump from −0.35/N to 0. Boltzmann results stay bit-identical through
    `models.fermi_dirac = false`, a separate path (6.2, Unit 14).

  - OLD (`inputs.cpp:83-88`, `device2d.py`, `device1d.cpp:303-316` [verified]): the contact neutrality used
    eg_kt = E_g/kT, the un-narrowed gap, while the bulk densities used ln_gn = ln(Nc/n_ie,eff). With band-gap
    narrowing on, the contact's minority density did not satisfy the bulk relation at the contact potential.
  - NEW: contact and bulk both use g_n and g_p.
  - REASON: one statistics relation. The two forms are the same without narrowing.

  - Not carried: incomplete ionization, which the legacy M13 phase 2 also had. It is independent of Fermi–Dirac,
    hydrogenic, and invalid above the Mott transition near 4e18 cm⁻³; it is deferred until requested (14.3). Also not
    carried: the legacy `ni_fd`, since ψ stays referenced to the Boltzmann intrinsic level, so n_ie keeps its meaning
    and the Unit 12 gate offset is unchanged. The legacy `band_diagram` is not carried either, as there is no
    band-diagram output.

**As built (Unit 15, materials and interface-local models; owner request, roadmap item 15):**
- `semiconductor.hpp` adds the other legacy material sets: `germanium_parameters`, `gallium_arsenide_parameters`,
  `indium_gallium_arsenide_parameters` (In0.53Ga0.47As), `silicon_carbide_4h_parameters`, and
  `algaas_parameters(x)` for Al_x Ga_1−x As, valid for 0 ≤ x ≤ 0.45 (`invalid_input` outside, where the gap is
  indirect). Every field a legacy set left at its dataclass default keeps the silicon value, as in the legacy:
  lifetimes and, where not given, Auger and the Canali exponents. All six sets validate; E_g, n_i and the depth below
  match 40-digit values at 300 K to 1e-12–1e-14. (Changed in the Unit 15 follow-up below: no Slotboom narrowing
  outside silicon, and AlGaAs built on GaAs.)
- `intrinsic_level_depth_eV(m, T)` = χ + (E_c − E_i), with E_c − E_i = E_g/2 + (kT/2) ln(Nc/Nv). Its step between two
  materials is the step of the potential at which each holds n = n_i, so it defines the band shift (6.3, Unit 15).
- `thermionic_emission.hpp`: `emission_velocity_cm_s(N, T)` = (kT/h)(2/N)^(1/3). This is the legacy
  sqrt(kT/(2π m)) with m recovered from N, simplified. Silicon's Nc gives the legacy 2.575e6 cm/s, matching 40
  digits to 1e-14. The legacy records it as about a factor 1.92 below the tabulated Richardson constant for silicon,
  because the density-of-states mass and the Richardson mass differ.
- Findings in the legacy sets, ported unchanged:
  - The 4H-SiC comment says E_g(300 K) ≈ 3.23 eV; its own Varshni numbers give 3.201 eV.
  - The AlGaAs comment says the conduction band takes about 85% of the gap step; its χ and E_g0 slopes give
    0.85/1.247 = 68%.

**As built (Unit 15 follow-up, materials; owner request "finish the follow-up" after the Unit 15 review):**
- `SemiconductorParameters` gains `radiative_cm3_s` (B), `ionization` (`IonizationParameters`: donor depth below E_c,
  acceptor height above E_v, degeneracies) and `richardson` (`RichardsonParameters`: effective Richardson constants, 0 =
  the density-of-states velocity). Validation: B, the levels and the Richardson constants ≥ 0, the degeneracies > 0.
- Values (not in the legacy, sourced): radiative B for Ge 6.4e-14, GaAs 7.2e-10, In0.53Ga0.47As 0.96e-10, 4H-SiC 1.5e-12
  (the bimolecular value), AlGaAs 1.8e-10 cm³/s (Ioffe; nextnano for GaAs); silicon 0 as in the legacy (its 1e-14 is
  negligible against SRH), so silicon results do not change. Dopant levels: silicon the legacy hydrogenic 45 meV (M13),
  GaAs Si donor 5.8 and C acceptor 26.3 meV, Ge P 12.0 and B 10.4 meV, 4H-SiC N 70 and Al 220 meV; InGaAs and AlGaAs 0
  (complete ionization; AlGaAs donors form DX centres, which one shallow level does not describe). Richardson constants
  0 everywhere (the legacy's choice); the legacy's Schottky table values are quoted in the header for users who set
  them.
- OLD / NEW / REASON, material defaults:
  - OLD: every non-silicon set carried silicon's Slotboom band-gap narrowing, and AlGaAs took silicon's defaults for
    every field it did not set (μ_min 92 and 47.7, silicon Auger and saturation).
  - NEW: no narrowing outside silicon (`no_bandgap_narrowing`), and AlGaAs is built on the GaAs set.
  - REASON: the Slotboom form and numbers are a silicon fit; a GaAs-alloy built on silicon's fields mixes two materials.
- `algaas_parameters(x, conduction_share)`: the conduction band takes `conduction_share` of the gap step (default
  0.85/1.247 = 0.68, the legacy χ slope; measured GaAs/AlGaAs offsets are nearer 0.62–0.65).
- `ionization.hpp`: `ionized_density(N, η, E/kT, g)` = N/(1 + g e^(η + E/kT)) and its η-derivative, overflow-free
  (legacy `ionized_eta_doping`); `ionized_neutral_equilibrium` solves n − p = N_D+ − N_A− by safeguarded Newton from the
  completely ionized root, under either statistics. Boron in silicon at 77 and 300 K and the 4H-SiC dopants at 300 K
  match 40-digit roots (legacy G7(b,c)).
- `radiative_recombination(n, p, E, B)` = B(np − E) with the partials of E; `emission_velocity_cm_s(m, carrier, T)` = A*
  T²/(q N) when the material sets A*, else the density-of-states velocity.

**As built (Unit 15b, insulators and interface traps; owner request, scope approved by the owner):**
- `physics/insulator.hpp`: `InsulatorParameters{eps_r}` and `Insulator::create` (eps_r finite and positive). Built-in
  set: SiO₂, 3.9 (the legacy `EPS_OX_R`); any other insulator is given by its permittivity (the owner chose this over
  sourced Si₃N₄, Al₂O₃ and HfO₂ sets).
- `physics/interface_traps.hpp` (not in the legacy, whose only trap model is `moscap`'s lumped q D_it, M14):
  - `TrapLevel` (donor- or acceptor-like, sheet density N_t, energy E_t − E_i, σ_n, σ_p), `TrapBand` (uniform D_it
    between two energies), `InterfaceTraps` (levels, bands, thermal velocities, default 1e7 cm/s);
    `check_interface_traps` against the semiconductor's gap at the device temperature.
  - Steady-state SRH occupancy f = (c_n n + c_p p1)/(c_n(n + n1) + c_p(p + p1)) and rate
    U/N_t = c_n c_p (np − n1 p1)/(…), c = σ v_th, with n1 = γ_n n_ie e^τ and p1 = γ_p n_ie e^−τ (γ the degeneracy factor
    of the node's own density, 1 under Boltzmann). At equilibrium n1/n = p/p1 = e^(τ − η), so f is exactly the Fermi
    function under either statistics (measured within 1.9e-15 over 1134 equilibrium states) and U vanishes.
    `trap_kinetics` returns f, 1 − f (computed directly, no cancellation) and U/N_t with exact partials, including the
    density dependence of n1 and p1.
  - Charge signs (owner's constraint): a donor-like trap is +q when empty, an acceptor-like trap −q when occupied.
  - Trap bands: composite 6-point Gauss–Legendre on panels no wider than kT (`trap_band_levels`). The occupancy's poles
    lie π kT off the real axis, so the rule converges fast: on a 1.08 eV band at 300 K (252 levels) the equilibrium
    occupied density is within 7.1e-16 of the band's density of the closed-form integral (40-digit references).
  - `fermi_occupancy(x)`: f and 1 − f from the exponential that cannot overflow, within 4 ulp of 40-digit values from
    x = −40 to 45.
  - The surface recombination velocity form is SRH with a mid-gap level per unit area (`srh_recombination` with
    τ = 1/s); a mid-gap trap with N σ v = s gives the same rate to 8 ε.

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
  - Heterojunctions (an edge between regions whose material parameters differ) were rejected as `invalid_input`
    until Unit 15 (6.3), which added band offsets and permittivity steps. Disconnected regions could already differ.
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
  charge-neutral potential and returns the result in V and cm⁻³ (6.2); since Unit 10 that is
  `results::EquilibriumResult` (6.6). The linear solver keeps its default configuration, including the
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
  - The convergence record (Unit 10: `results::ConvergenceRecord`) keeps the smallest pivot ratio.
- `DriftDiffusionModels{doping_mobility, srh}` (Unit 11: `PhysicsModels`, with `auger` and `bgn`) are the legacy `Models`
  flags of the baseline diode. OLD / NEW / REASON:
  with `doping_mobility` off, the legacy used `mu_max` at 300 K whatever the temperature; NEW uses Caughey–Thomas at N = 0,
  so μ_max(T); REASON: the same model at all temperatures (identical at 300 K).
- `solve_bias(device, bias_V, options, initial)` starts from the given state (a sweep's previous point) or from the
  equilibrium solution, sets the contact nodes to their biased values (legacy) and returns the potential and densities,
  terminal currents, edge currents and the convergence history; since Unit 10 that is `results::BiasPoint` (6.6), and
  `solve_bias` is a sweep of one point (6.9).
- Not carried from the legacy solver: the line search and generation-strength stages (impact ionization and BTBT only),
  Robin contacts with surface recombination velocity, the energy-balance block, and lagged field mobility. (The ln(n_ie)
  term in δ was added with band-gap narrowing at Unit 11, section 5; field mobility, not lagged, at Unit 13, below.)
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

**As built (Unit 13, field-dependent mobility; owner request, roadmap item 13):**
- `SemiconductorParameters` gains `electron_saturation` and `hole_saturation`, `CanaliParameters{v_sat_cm_s, beta}`,
  with the legacy silicon values (1.07e7 cm/s, β = 2; 8.37e6 cm/s, β = 1), temperature independent as in the legacy.
  Validation: v_sat > 0 and β ≥ 1.
  - OLD: the legacy required β > 0.
  - NEW: β ≥ 1.
  - REASON: below 1 the derivative dμ/dE is infinite at E = 0, and the Jacobian needs it. The legacy's own values are 1
    and 2.
- `physics/field_mobility.hpp`: `canali_mobility(mu0, E, CanaliParameters)` returns the legacy μ0/[1 +
  (μ0E/v_sat)^β]^(1/β) and its exact derivative in E (the right-hand one at E = 0: 0 for β > 1, −μ0²/v_sat for β = 1).
- `PhysicsModels::field_mobility` (off by default, as in the legacy), in the drift-diffusion assembler only. The run
  record lists it (not in the quasi-static sweep, which has no transport).
- OLD / NEW / REASON, the field-mobility coupling:
  - OLD (legacy `Device1D`, `field_mobility`): per edge |E| from ψ, averaged to the nodes, Canali applied per node, then
    the harmonic mean on the edges. The mobility was lagged: recomputed from the iterate before each Newton step, with
    no Jacobian term. Legacy `Device2D/3D` refused the model.
  - NEW: an edge-local model (section 5, item 4). Each edge's low-field mobility (the harmonic mean, as before) is
    scaled by the Canali factor of the edge's own field E = V_T|ψ_b − ψ_a|/length, the field component along the edge
    and hence parallel to the edge's current. dE/dψ enters the Jacobian, which stays inside the edge pattern, so Newton
    stays quadratic. It works unchanged on the mesh graph in 1D, 2D and 3D.
  - REASON: a lagged mobility costs the quadratic convergence and is not checked by the FD-Jacobian gate. The node
    average mixes the fields of edges in different directions on a 2D/3D graph. On a uniform field the two forms agree.
  - The driving field stays the electrostatic one, as in the legacy. Using the quasi-Fermi gradient instead was tried
    and rejected: with it, Newton on a 1e19/1e17 diode failed at the first forward and reverse step, where the
    electrostatic form converges. The minority quasi-Fermi potentials (ln of densities near 1e-14) make the edge factors
    erratic during the iteration.
- **Finding, damping with field mobility in 2D:** on the legacy MOSFET at V_G = 1 V, the first 50 mV drain step diverges
  under the legacy cap of 5 V_T on the ψ correction (the corrections pass 1e15 within 13 iterations). It converges under
  1 V_T, in 10–24 iterations per point of a 0.1 V sweep. Gate ramps at V_D = 0 converge with the default cap. The
  mechanism was not isolated. The defaults are not changed (they are the legacy ones); a line search or a field-aware
  damping are candidate fixes, when a unit asks for one.
- Measured effect: J(0.5 V) of the legacy 1e17 diode is 0.77% lower with the model on (the Unit 9 gate, with the model
  off, is unchanged). The MOSFET drain current at V_G = 1 V falls to 0.968 of the constant-mobility value at V_D = 0.1 V
  and to 0.735 at 1 V.

**As built (Unit 14, Fermi–Dirac statistics in the assemblers; section 5 for the physics):**
- `PhysicsModels::fermi_dirac` is off by default, like the legacy `fd`. Both assemblers read it, and the run record
  lists it for both equation sets.
- `EquilibriumPoisson` takes n, p and dn/dψ from the selected statistics; its contact values and charge-neutral guess
  are the Fermi–Dirac root. It is the legacy `solve_equilibrium_fd` without the η clamp.
- `DriftDiffusion` uses the legacy ν-factor scheme (M13 plan 3.2bis). δ_n gains ln γ_n,b − ln γ_n,a and δ_p loses
  ln γ_p,b − ln γ_p,a, each node's γ taken from its own density by `fermi_dirac_degeneracy`:
  - Then n_b/n_a = e^δn at equilibrium, so the flux vanishes there.
  - In the continuum limit the flux is μ n ∇φ_n: the generalized Einstein relation is implicit.
  - The Jacobian's density columns gain d flux/dδ times d ln γ/d density.
  - SRH and Auger drive n p towards n_ie² γ_n γ_p, with partials. Auger takes physical densities, so E = Ns² E′
    and dE/dn = Ns dE′/dn′.
  - `state_from_potential`, the contacts and the gate rows go through the same statistics. The gate's intrinsic-level
    offset (6.4) is unchanged, because the gauge is.
- **Boltzmann path bit-identical:** the same probe was built against `main` and against this unit (a scratch test, not
  committed). It hashes every potential, density, current and Newton record of an equilibrium solve and a nine-point
  sweep (0–0.8 V and −2 V) for four configurations: 1e17/1e17 with default models; 1e19/1e18; 1e18/1e16 with field
  mobility; and 1e17/1e17 with every model off. All four hashes are equal.
- Measured gates:
  - FD-Jacobian:
    - Poisson: 7.6e-10 / 6.9e-9 / 9.9e-9 in 1D / 2D / 3D on a 1e20/1e17 diode.
    - Drift-diffusion with every model on: 2.5e-7 / 2.2e-6 / 5.3e-6.
    - The degeneracy part alone (J(FD) − J(Boltzmann) against FD of the residual difference, majority columns):
      1.5e-7.
    - The Fermi–Dirac recombination part alone, at a near-equilibrium uniform 1e20 state where E's partials are
      comparable to p: 2.2e-11.
  - Detailed balance: every edge of the 1e20/1e17 junction, with band-gap narrowing, is within 1.3% of its rounding
    bound at the equilibrium state.
  - Legacy G7(a), uniform 1e20 in 1D and 2D: n = N_D and p = 0.3496, matching the 40-digit root to 1e-11. That is a
    third of the Boltzmann minority density.
  - Legacy G4(d), 1e20/1e17 without narrowing: V_bi = 1.03688 V, the root pair to 1e-13, and 28.3 mV above
    Boltzmann.
  - Legacy G4(c), generalized mass action at every node: 6e-15.
  - Legacy G6(b), 1e16/1e16: densities deviate by 1.24e-4 = δ (gate 3δ); currents at 0.5 V by 1.1e-4 (gate 20δ).
  - p+ 1e20 / n 1e17 forward sweep, all models: the electron current injected into the degenerate p+ side falls by its
    hole degeneracy factor γ_p = 0.3350 (measured 0.33538–0.33545 at 0.2–0.7 V). The hole current into the n side is
    unchanged to 7e-4, within that side's own δ = 1.2e-3. The diode current falls 1.7–2.4%; Newton takes 5–6
    iterations per point.
  - Legacy G7(d), the quasi-static MOS-C at N_A = 1e18 and 5 nm: C_max falls from 0.981 C_ox to 0.960 C_ox, −2.1%
    (legacy band 2–30%, so near its lower edge).
- **Finding, the zero-bias drift-diffusion step:** from the Fermi–Dirac equilibrium, the zero-bias solve takes one
  Newton step of 1e-10 (relative) and stops.
  - Cause: the equilibrium potential leaves a few ulp of charge per node, since ψ ≈ 22 resolves e^ψ only to 22 ε.
  - The minority electrons of the 1e20 region absorb that step at up to 2e-11 relative. That is solver accuracy, not
    the scheme, which is why detailed balance is gated at the equilibrium state itself.
- **Finding, test design:**
  - Relative-step differences cannot resolve the minority-density columns: rows that also carry majority fluxes of
    order 1 change by less than their rounding (the Unit 13 finding). The degeneracy-part check therefore uses the
    majority columns, and the full gate covers the rest with its absolute step.
  - Neither the full gate nor the degeneracy-part check sees the equilibrium product's partials (the Unit 9 finding),
    hence the separate near-equilibrium recombination check.
  - A mutation that dropped those partials survived until that check was added; it also exposed a physics test that
    skipped its comparison based on the analytic partial under test.
- Mutation checks: each of the nine below fails at least one test.
  - Dropping the electron density chain of the degeneracy term.
  - Dropping ln γ from δ_n.
  - The wrong sign of the hole degeneracy term.
  - A Boltzmann charge derivative in the Fermi–Dirac Poisson Jacobian.
  - The equilibrium product without its partials.
  - The degeneracy derivative without its 1 + d ln γ factor.
  - Boltzmann contacts under Fermi–Dirac.
  - The table's curvature dropped.
  - Sommerfeld cut to two terms.

**As built (Unit 15, heterojunctions in the assemblers; 6.3 for the model):**
- `make_scaled_device` no longer rejects heterojunctions. It adds the band shift per node, the harmonic-mean edge
  permittivity, an interface flag per edge, and the shift in each gate node's offset.
- `EquilibriumPoisson` and `DriftDiffusion` evaluate the statistics at ψ + s, and set contacts and the neutral guess
  to η₀ − s. The SG driving terms are δ_n = Δψ + Δs + Δln n_ie (+ Δln γ_n) and δ_p = Δψ + Δs − Δln n_ie
  (− Δln γ_p).
- `PhysicsModels::thermionic_emission` is off by default and applies to drift-diffusion only. On an interface edge
  the fluxes are those of `thermionic_flux.hpp` (legacy M33-S2):
  - Jn = K(n₂g₂ − n₁g₁), with g₁ = min(1, e^u), g₂ = (Nc₁/Nc₂) min(1, e^−u), u = δ_n − ln(Nc₂/Nc₁); holes mirror
    this.
  - K is the harmonic mean of the ends' emission velocities times the scaled interface area.
  - Detailed balance holds exactly for either step sign, and with Fermi–Dirac too, because δ carries Δln γ.
  - The run record lists the switch for drift-diffusion.
- `DriftDiffusion::create` refuses field mobility on a device with an interface edge.
- **Newton measure floor (changes the Unit 9 criterion):**
  - OLD (Unit 9): max(|dψ|, |dn|/n, |dp|/p).
  - NEW: each density is measured against max(n, 10⁻²⁰ × the largest density of the state).
  - REASON, measured:
    - A GaAs/AlGaAs junction at 300 K has minority electrons of 2e-11 cm⁻³ in the AlGaAs, 2e-30 of the largest
      density. Newton's residual converged (1e-13) while |dn|/n wandered between 0.01 and 0.1 for 100 iterations,
      with every model combination, plain Boltzmann included.
    - The linear solve resolves a density correction only to about ε max|dx|, and near convergence max|dx| ≈
      ε n_max. So densities below about 5e-24 n_max cannot reach 1e-8 relative to themselves.
    - With the floor, all five combinations converge through 0.6 V forward in 1–21 iterations per point.
    - Silicon at 300 K stays near 1e-14 of its largest density, above the floor, so its iterates are unchanged
      (6.3, bit-identity).
- Measured gates:
  - FD-Jacobian on a Si / Si(χ + 0.3, ε_r + 1.3, E_g0 + 0.1) p-n junction, every column, in 1D / 2D / 3D: 1.2e-9 /
    3.3e-9 / 1.8e-9, also with thermionic emission, and 1.9e-9 / 2.5e-9 / 4.9e-9 with Fermi–Dirac as well.
  - The same junction under Poisson alone: 1.9e-9.
  - The thermionic part alone (J(TE) − J(SG) against FD of the residual difference): 2.0e-9.
  - The thermionic edge current equals hmean(v) L_D/D0 × (g₂ − g₁) to 1e-13.
  - The 1D Poisson row across a χ and ε step, written out, to 1e-12.
  - Contacts and the neutral guess move by exactly the step, while the carriers at the guess equal the
    homojunction's.
- **Thermionic emission, measured:**
  - On the legacy isotype fixture (n-Si / n-Si(χ), 1e17, 0.1 V), J(TE)/J(DD) is 1, 0.192 and 0.034 for steps of 0,
    0.15 and 0.30 eV. The legacy S2 gates hold: emission limits the current, and the limit deepens with the step.
  - At no step the two materials are one, so there is no interface and the ratio is exactly 1.
- Mutation checks: each of the ten below fails at least one test.
  - The band-shift sign flipped in Poisson.
  - Drift-diffusion contacts without −s.
  - The SG driving term without Δs.
  - An arithmetic instead of harmonic edge permittivity.
  - The thermionic g₂ without its Nc ratio.
  - The thermionic slope of g₁ dropped.
  - The gate offset without s.
  - No density floor in the Newton measure.
  - The emission velocity with (1/N)^(1/3).
  - The field-mobility refusal removed.

**As built (Unit 15b, insulators in the assemblers; 6.4 for the device data, section 5 for the models):**
- An insulator node has no carriers. In `DriftDiffusion` its continuity rows are F = n and F = p (identity rows, n = p
  = 0 held exactly, owner's constraint) and its Poisson row the box row with no charge; in `EquilibriumPoisson` the box
  row with no charge. An edge with an insulator end carries Poisson flux only (`ScaledEdge::carriers`), so no carrier
  crosses a semiconductor-insulator edge: its edge currents are 0 and no semiconductor row reads an insulator density.
  The edge permittivity is the harmonic mean, as at a heterointerface; the interface lies at the edge's midpoint and
  each control volume is of one material (`mesh::straddle_interface`), so no mesh change was needed.
- An electrode node's Poisson row is Dirichlet, ψ − ψ_E, ψ_E = V/V_T + (depth of the reference material − φ_m)/V_T
  (the vacuum level is −ψ V_T + depth). It is the lumped gate's ψ_G with the oxide meshed. The electrode's charge is
  the displacement flux leaving its nodes along the edges to other nodes, reported in `gate_charge` with the gates'.
- `assemble/interface_edges.hpp` (`InterfaceEdges`, shared by both assemblers as `GateNodes`): each mesh edge joining
  the insulator and semiconductor regions of a declared interface with charge, traps or recombination is two half-edges
  in series (g_i = 2 et_i c, g_s = 2 et_s c), and the interface potential ψ_I between them obeys Gauss's law on the
  edge's interface patch (area A, a = A/L_D^(D−1)): g_i(ψ_I − ψ_i) + g_s(ψ_I − ψ_s) = a/(Ns L_D) (Q_f + Σ N_k q_k). The
  insulator row takes g_i(ψ_I − ψ_i) and the semiconductor row g_s(ψ_I − ψ_s) in place of the edge's flux (with no
  charge this is the harmonic-mean flux); the semiconductor's continuity rows take ∓ a Ns/(R0 L_D)(Σ N_k r_k + r_s). In
  `EquilibriumPoisson` the occupancy is the Fermi function of τ − (ψ_I + s); in `DriftDiffusion` the carriers keep their
  quasi-Fermi levels across the semiconductor half-cell (Boltzmann: n_I = n_s e^(ψ_I − ψ_s), p_I = p_s e^(ψ_s − ψ_I);
  Fermi–Dirac through the node's reduced energy), and the occupancy is the SRH one at (n_I, p_I). Q falls as ψ_I rises
  and lies between the fixed charge with every trap charged one way or the other, so ψ_I is the one root of a monotone
  equation in a known bracket, found by safeguarded Newton to rounding on each evaluation. ψ_I is eliminated (the global
  unknowns are unchanged); its partials in the edge's end unknowns follow from the local equation (implicit function),
  so the Jacobian is exact. In drift-diffusion the insulator node's Poisson row then reads the semiconductor node's n
  and p: the pattern gains those two entries on interface edges only. A Dirichlet end (electrode, ohmic contact) takes
  no terms; an electrode at an interface edge's end has the half-edge flux as its charge.
- OLD / NEW / REASON, where the interface terms act (within Unit 15b, owner's request after the first commit `daf4de7`):
  - OLD: the fixed charge and traps acted at the first semiconductor node, h/2 from the interface, at its own n, p and ψ.
  - NEW: at the interface potential ψ_I and the interface densities, per interface edge.
  - REASON: the node placement was first order in h, about q D_it E_s h/2 in the trapped charge and (h/2)/ε_si in the
    fixed charge's shift, while the rest of the MOS solution is second order. Measured with D_it = 2e12 against the
    exact solution: 3.8e-4 V (1200 nodes) halving with h, now 1.77e-5 V falling by 4.04–4.08 per halving.
- Newton: an insulator node's densities are neither measured nor updated (they stay 0).
- The reference material (scaling, band shift) is the lowest semiconductor node's (`Device::reference_node`): node 0
  when there is no insulator, so nothing changes for existing devices.
- Band diagram: NaN on insulator nodes (the model has no insulator bands).
- Gates (all measured):
  - FD Jacobian in 1D, 2D and 3D, Boltzmann and Fermi–Dirac, electrode biased, fixed charge, two trap levels, a trap
    band and surface recombination: at most 1.7e-7 (Poisson) and 1.2e-7 (drift-diffusion); the interface part on its
    own (difference of the Jacobians with and without the terms) at most 2.0e-8.
  - On an equilibrium state the drift-diffusion interface terms are the equilibrium ones: trapped charges equal, Poisson
    rows equal to rounding, continuity rows at most 1.9e-14, in 1D/2D/3D under either statistics.
  - A mid-gap trap with N σ v = s (N small, so its charge does not move ψ_I) gives the velocity form's continuity rows
    exactly (0 difference measured).
  - Meshed MOS-C (legacy fixture: N_A = 1e17, 5 nm SiO₂, n+ poly) against the exact first-integral solution, gate charge
    over −2 to 2 V: second order, 5.0e-4, 1.24e-4, 3.07e-5, 7.6e-6 V for 300 to 2400 silicon nodes; the lumped oxide
    6.6e-4 to 1.08e-5; meshed against lumped 1.16e-3 to 1.8e-5.
  - Fixed charge: the flat-band shift equals −q Q_f/C_ox to 6e-13 (Q_f = 5e11 and −1e12 cm⁻²); the curve's error
    2.4e-5 and 3.5e-5 V (it was 9.4e-5 and 2.1e-4 V with Q_f at the node).
  - Uniform D_it = 2e12 cm⁻² eV⁻¹ neutral at flat band: against the exact solution 2.94e-4, 7.20e-5, 1.77e-5, 4.39e-6 V
    for 300 to 2400 nodes (second order, ratio 4.04–4.08); the trapped charge within 1.1e-12 C/cm² of the exact Fermi
    integral; the legacy M14 stretch-out (the ported legacy solve with its D_it term) within 0.17 mV of surface
    potential where the surface Fermi level is mid-gap, against 236 mV without traps.
  - Gauss's law with electrode, fixed charge and traps (bound 1e-9 relative), quasi-static and drift-diffusion, also with
    a one-cell oxide whose electrode is the interface edge's end; y-uniform 2D and 3D reproduce 1D (bound 1e-10);
    drift-diffusion equals the quasi-static state below threshold; a metal electrode moves the curve rigidly by φ_m − χ.
  - Surface recombination: a 2D p+n diode under oxide with s = 0, 1e3 and 3e3 cm/s against the analytic current (lowest
    transverse mode, k tan(k w) = s/D) within 8.3e-4 (2.7e-4 at s = 0).
  - MOSFET (legacy fixture, `[.mosfet]`): meshed oxide against lumped, thresholds 0.21 mV apart, currents from V_G =
    −0.4 V up within 0.58%.
  - Without interface charge, traps or recombination nothing changed: devices without insulators hash identically to
    `main` (seven runs), meshed-oxide devices without interface terms identically to `daf4de7` (1D and 2D, quasi-static
    and drift-diffusion).
  - Cost (Release, best of three; run-to-run noise on this machine is 10–20%): with traps, a 1D quasi-static sweep of 81
    points 0.067 → 0.077–0.083 s, 1D drift-diffusion 16 points 0.151 → 0.159–0.180 s, 2D (8 columns) quasi-static 21
    points 0.86 → 0.88–1.03 s, 2D drift-diffusion 7 points 3.81 → 4.49–4.66 s; Newton iterations 350 → 350, 116 → 114,
    121 → 121, 65 → 68. The sparse factorization still dominates.
- Known limits:
  - Trap occupancy is steady state in a bias solve; the transient run follows its dynamics (Unit 21, below).
  - A coarse semiconductor half-cell at the interface makes the interface densities steep in ψ_I (n_I = n_s e^(ψ_I − ψ_s)):
    with a one-cell 5 nm oxide (a 5 nm half-cell) in accumulation, the drift-diffusion Jacobian's pivot ratio falls to
    1.2e-10 (6.9e-5 with the terms at the node), near the default singularity check of 1e-11; on the legacy mesh
    (0.0125 nm half-cell) it is 1.5e-5. Accuracy needs a fine interface mesh anyway.
  - Not in 15b (owner's exclusions): tunnelling through insulators, poly-gate depletion, bare-surface recombination
    velocity, per-cell regions in the mesh.
- Mutation checks, each caught (26, after the interface-potential change): the donor charge as occupied, the trap
  charge's n partial, the trap rate's p partial, p for p1 in the occupancy, the half-edge conductance without its
  factor 2, the interface area halved, carrier flux across semiconductor-insulator edges, the electrode potential's
  sign, the electrode charge's sign, the charge of an electrode at an interface edge's end (equilibrium and
  drift-diffusion), the fixed charge missing in drift-diffusion, the band quadrature weight, semiconductor parts joined
  through oxide in the topology check, insulator doping unchecked, the recombination weight without Ns, the thermal
  velocity ignored (missed at first: the trap-against-velocity test compared rows on their own scale, and now compares
  the recombination term on its), the run record without the electron thermal velocity (missed in the first commit; the
  run-identity test was extended), a charge on insulator rows, the occupancy's sign of τ − η, ψ_I's implicit
  derivative dropped, its denominator without dQ/dψ_I, the interface electron density shifted the wrong way, the local
  solve stopped after one step, an interface edge also taking its own flux, the insulator row without its n column.

**As built (Unit 21, transient drift-diffusion; owner request, scope approved by the owner; legacy `transient.py`,
`transient2d.py`, `transient3d.py`):**
- **Time steps in the assembler** (`DriftDiffusion::TimeStep`, `evaluate`/`residual` with a step). A BDF step from
  earlier states to t_new solves the steady rows plus a storage term on each semiconductor node off the ohmic contacts:
  electrons − V_i r (S_n − c_n), holes + V_i r (S_p − c_p), r = 1/(β h) with h in units of t₀ = Ns/R0 = L_D²/D0 (the
  time scale of the scaled rates, the legacy `_time_scale`), c the step's combination of the earlier storages. Poisson,
  contact, electrode and insulator rows have no storage. Every steady model carries over: SG and thermionic fluxes,
  field mobility, Fermi–Dirac, band-gap narrowing, lumped gates, meshed insulators.
- **What a node stores.** S_n = n and S_p = p; with incomplete ionization S_n = n − N_D⁺ and S_p = p − N_A⁻, the
  carriers bound to the dopants counted with the free ones. Ionization is instantaneous in the model, so storing n
  alone would create or destroy charge whenever N_D⁺ moves; with S the charge balance is exact (measured: the total
  currents sum to 1e-11 of the largest with the model on, and the model changes the current).
- **Interface trap dynamics** (the 15b exclusion). Per level, df/dt = occ − f D (occ = c_n n_I + c_p p1,
  D = c_n(n_I + n1) + c_p(p_I + p1)); a BDF step f − c = β h (occ − f D) is linear in f, so
  f = (c + k occ)/(1 + k D), k = β h Ns, is eliminated in closed form inside the per-edge ψ_I solve: the global
  unknowns stay as they are and the Jacobian stays exact. f is a weighted mean of c and the steady occupancy occ/D, so
  it lies in [min(c, 0), max(c, 1)] and the trapped charge within the bounds of those, which keep the local solve
  bracketed; as k grows it is the 15b steady state. The electron row takes the net electron capture
  c_n(n_I(1 − f) − n1 f) and the hole row the net hole capture c_p(p_I f − p1(1 − f)); they differ by N df/dt. The
  history of each (interface edge, level) is a trap slot (`InterfaceEdges::slot_offset`).
- **Contact charges** (`DriftDiffusion::contact_charges`): the displacement flux leaving each contact's nodes along the
  edges to nodes outside it (an electrode as 15b, a gate as 12, now also an ohmic contact). Gauss's law makes the
  contact charges balance the charge inside the device (carriers, dopants, fixed and trapped interface charge); the
  residuals of the Poisson rows are exactly that balance (assemble test, 1e-12).
- **Conduction current** of a contact (`conduction_currents`): the steady `terminal_currents` plus, on an interface edge
  whose semiconductor node is on an ohmic contact, the current that charges its traps from the contact (rate_p − rate).
- **Integrators** (`solve/transient.hpp`): backward Euler (the legacy's) and variable-step BDF2, the default:
  y_{n+1} − ((1+ω)² y_n − ω² y_{n−1})/(1+2ω) = h (1+ω)/(1+2ω) F(y_{n+1}), ω = h_n/h_{n−1}; second order, L-stable, one
  Newton solve per step. The run starts, and restarts after every waveform breakpoint, with two backward-Euler steps (the
  second with an error estimate), then BDF2. TR-BDF2 was considered and not used: its trapezoidal stage has an explicit
  half that rings on stiff modes (dielectric relaxation is about 1e-13 s; fast traps likewise), which can take densities
  negative or occupancies far outside [0, 1] within a step.
- **Step control.** Milne's device: the local error is C/(P − C) times the difference between the solution and the
  polynomial through the earlier states (linear for backward Euler, C = −h²/2; quadratic for BDF2,
  C = −(1+ω)² h³/(6ω(1+2ω))), measured as |e_ψ| (units of kT/q), |e_n|/(n + n_ref), |e_p|/(p + p_ref) and |e_f| on trap
  occupancies, against rtol (default 1e-3; n_ref 1e10 cm⁻³). A step whose estimate exceeds rtol is retried at
  h (0.9/ratio)^(1/(order+1)) (at least 0.2 h); an accepted one sets the next h by the same rule within [0.2 h, 2 h]
  (2 is below BDF2's zero-stability limit 1 + √2). A step whose Newton solve fails is retried at h/4. Below dt_min the
  run stops with non_convergence and keeps its steps. Steps land exactly on waveform breakpoints, output times and
  t_end; the first step after a start or a breakpoint is at most dt_initial and has no estimate. Fixed steps
  (`adaptive = false`) take dt_initial throughout.
  - OLD / NEW / REASON: OLD, the legacy grew the step by 1.5 after a Newton solve of few iterations and halved it on
    failure; NEW, error-controlled steps; REASON, the number of Newton iterations says nothing about the time error.
- **Terminal currents.** Total = conduction + displacement, the displacement current being the step's own BDF
  difference of the contact charge, (Q_c − c_Q)/(β h). Every quantity stored (carriers, bound carriers, trapped and
  contact charge) is differenced by the same formula, so the total currents of all contacts sum to zero to the solve's
  tolerance, and under backward Euler Σ h I_disp telescopes to the change of each contact's charge.
  - OLD / NEW / REASON: OLD, the legacy reported Jn + Jp on the contact edge only; NEW, conduction plus displacement
    current; REASON, the conduction current alone is not conserved in a transient (a gate carries only displacement
    current).
- **Waveforms** (`solve/waveform.hpp`): piecewise linear (two corners at one time are a jump; the waveform is
  right-continuous, and a step ending on the jump uses the left value), with the legacy step, ramp and pulse built on
  it (a pulse may have rise and fall times), and sine (SPICE SIN, with delay and phase). Corners and the sine's delay
  are breakpoints. Every contact, gate and electrode included, takes one.
- **Starting state:** the steady state at the waveforms' values at t = 0 (left side), from thermal equilibrium or a
  given guess; the traps start at their steady occupancy, so the initial state is consistent.
- Gates (all measured):
  - FD Jacobian of the time-step system (steps of 1 ps and 1 ns, histories off the state, trap histories in
    [−0.05, 1.05]), 1D/2D/3D, Boltzmann and Fermi–Dirac, with incomplete ionization: at most 1.7e-7; the trap part on
    its own at most 3.5e-7 (difference step 1e-6: at 1 ps the house step's rounding shows, 3e-6, falling as the step
    grows, so it is rounding, not the Jacobian).
  - A step of 1e12 s gives the steady rows and occupancies (difference 0 and 2.2e-16); the capture imbalance of the two
    rows equals r ΔQ_trap to 1e-6 (rounding).
  - Constant bias from steady state: the current stays the steady one to 1.9e-15, no displacement current, 22 steps
    from 1 ps to 1 µs, each at most twice the one before. A diode stepped to 0.5 V ends at sweep_bias's state
    (bounds: current 1e-8, potential 1e-9 V, densities 1e-7; the current agrees to the eleven digits printed).
  - Conservation: a 2D MOS structure with an electrode, two ohmic contacts (one on an interface edge), fixed charge,
    a trap level, a trap band and surface recombination under a gate ramp and a pulse on a contact: the total currents
    sum to 3.1e-10 (backward Euler) and 2.3e-10 (BDF2) of the largest current; under backward Euler Σ h I_disp equals
    each contact's charge change to rounding. The diode turn-on: 4.1e-8. Traps on an ohmic node (a 1D oxide whose one
    silicon node is the contact): under backward Euler the contact's conduction current integrates to the trapped
    charge's change within 1e-9, and the currents sum to 3.8e-12.
  - Order (fixed steps, turn-on ramp, against BDF2 at 5120 steps): backward Euler error ratios 1.92–1.99 per halving,
    BDF2 3.33, 3.60, 3.78, 3.90 (its two backward-Euler starting steps fade from the ratio).
  - Error control: over a trapezoid with a jump, the charge error falls from 4.1e-3 to 9.9e-4 between rtol 1e-3 and
    1e-4 and the step count grows 2.0 times (10^(1/3) = 2.15); the jump's two sides are steps' ends and the integrator
    restarts there.
  - MOS-C RC response (10 nm oxide, 100 µm of 1e16 p-silicon, accumulated, 1 mV step): the gate charge relaxes with
    τ = R(C + C_g), R = L/(q μ_p N_A), C the quasi-static capacitance, C_g = ε_si/L: 0.092% (meshed oxide) and 0.093%
    (lumped gate) short, the local rate constant to 1e-4 over 0.05–3 τ.
  - Dielectric relaxation (100 µm oxide over 10 µm of 1e16 silicon, τ 1.03 times ε/(q μ_p N_A) = 1.59 ps): 0.34% short.
  - Trap emission: an acceptor level 0.3 eV below midgap (1e9 cm⁻²) after a gate step to weak inversion fills at
    B = c_p(p_I + p1) + c_n(n_I + n1) within 3.2e-4; its charge ends at the steady value to 2.2e-8.
  - The legacy diode turn-off (legacy `transient.py` built from the reference checkout, its p+n fixture, fixed steps of
    t_t/100): the anode conduction current equals the legacy's to 1e-13 on most steps and 3.7e-5 at worst (where it has
    fallen to 1.7e-4 of the forward current; both Newtons stop at 1e-8 of the densities, at different iterates). Legacy
    finding: from step 184 the legacy current repeats −1.4174e-6 A/cm² to the end, though the diode relaxes to
    equilibrium at 0 V; near equilibrium its merit line search finds no decrease at rounding level, takes no step
    (λ = 0) and reports the unchanged state as converged. NiTCAD keeps decaying (1.6e-14 A/cm² at 3 t_t).
  - y-uniform 2D and 3D reproduce the 1D transient to 2.3e-13; the meshed-oxide MOSFET (`[.mosfet]`, 45 × 25) after a
    gate step 0 → 1 V at V_DS = 0.1 V reaches the Id–Vg sweep's drain current within 1e-6, currents summing to 1.3e-12.
  - Steady paths unchanged: devices without insulators, meshed oxides without interface terms and meshed oxides with
    fixed charge and a D_it band (1D and 2D, quasi-static and drift-diffusion) hash identically to `main`, run digests
    included.
  - Cost: the steady paths are unchanged (a null-pointer test per node). A step costs one Newton solve under backward
    Euler and BDF2 alike. The 1D runs of the tests take 22 to 2948 steps; the 2D MOS structure 795 (backward Euler)
    and 331 (BDF2) steps, 4.5 s for both; the 45 × 25 MOSFET gate step 146 steps, about 50 s with its Id–Vg sweep
    (Release). The five slowest transient tests (81–280 s each in Debug, 16 s together in Release) carry the hidden
    tag `[.transient]` and run in Release only (ctest `solve_transient`), as the MOSFET gates do.
- Known limits:
  - The first step after a start or breakpoint has no error estimate; it is taken at dt_initial (default t_end/1e6),
    which the caller should keep small after a jump.
  - BDF2's history combination may leave a trap occupancy slightly outside [0, 1] (2.1e-8 measured); the bounds of the
    local solve allow for it.
  - Not in Unit 21 (owner's exclusions): AC small-signal (22), circuits and mixed mode, optical or external generation,
    dopant ionization dynamics, bulk trap dynamics, thermal transients (23), tunnelling, transient runs of the
    quasi-static equations.
- Mutation checks, each caught (27): the storage term's sign, its Jacobian, the bound donor electrons not stored, the
  trap history ignored, the occupancy's n partial, the hole capture without p1, the hole row taking the electron
  capture, the trap bounds without the history (missed at first: a random history never takes the whole interface
  charge out of the [0, 1] bounds; a test with every history at 1.2 and −0.2 was added), the trap weight without Ns,
  BDF2's β, its history sign, BDF2 never used, the displacement current's sign, the displacement current differenced by
  backward Euler under BDF2, the charge sign of an ohmic node on an interface edge, the ohmic contacts' charge left out,
  the traps' charging current at an ohmic node dropped (missed at first: in the 2D structure it is below the sum's
  bound; the ohmic-node test was added), backward Euler's error constant, BDF2's predictor on two points, the density
  error not relative, no restart at breakpoints, no landing on output times, the steady trap occupancy stored for the
  step's, the right side of a jump as a step's end bias, the waveform's exact left value at a jump (missed at first:
  the test's corner values interpolate exactly; one that rounds was added), the run record without the waveforms, the
  step growth unlimited (missed at first: the constant-bias test bounded the step count from above only; it now checks
  each step's growth).

### 6.3 Device description

A `device` is plain data: regions (geometry in the mesh's coordinates), doping per node or region,
material per region, and contacts (6.4). It contains no solver state and can be built without a solver
present. Per-region materials are in scope from the start of the data model; heterojunction physics came
with Unit 15 (below).

**As built (Unit 15, heterojunctions; owner request, scope chosen by the owner: semiconductor heterojunctions with
thermionic emission, meshed insulators in a later unit):**
- An edge whose two ends are in materials with different parameters is a heterointerface. The interface lies at the
  edge's midpoint, so each node's control volume is of its own material (the legacy's node-material model). The
  edge's permittivity is the harmonic mean of its ends', which is the series permittivity of the two half-edges, so
  D is continuous (legacy `_eps_tilde_edge`).
- **Band shift.** Each node carries s = (depth − depth₀)/V_T, with depth the `intrinsic_level_depth_eV` of the
  node's material and depth₀ that of node 0's material.
  - The carriers see ψ + s: n = n_ie e^(ψ+s) and p = n_ie e^−(ψ+s) under Boltzmann statistics, and η = ±(ψ + s)
    under Fermi–Dirac.
  - Contacts and the neutral guess are η₀ − s.
  - The SG driving terms gain s_b − s_a.
  - A gate's electrode potential becomes ψ_G − s at its nodes.
  - ψ is then the electrostatic (vacuum-level) potential, continuous across the interface, and s carries the band
    offsets.
  - On a device of one material s ≡ 0, and every expression reduces to the previous one operation for operation.
- The legacy M33 gates were ported, plus an analytic one. All measured:
  - **Exact interface:** the first integral of Poisson's equation on each side of an abrupt Boltzmann
    heterojunction, D² = 2qV_T ε (n_b(e^u − 1 − u) + p_b(e^−u − 1 + u)), with D continuous and the bulk potentials
    from a common Fermi level. For n-GaAs/N-Al0.3Ga0.7As, p-GaAs/N-AlGaAs (V_bi 1.524 V), N-AlGaAs/n-GaAs and a
    silicon χ+ε step:
    - V_bi matches to 1e-12.
    - The interface potential is within 1.4e-5 V.
    - D is within 7.6e-4.
  - **Legacy G1, detailed balance per carrier:** for χ steps of −0.3 to +0.3 eV and gap steps of ±0.2 eV, each with
    thermionic emission off and on, every edge at the equilibrium state is within 64 ε of its one-sided terms, and
    the legacy's 1e-7 A/cm² floor holds on the published zero-bias currents.
  - **Legacy G2:** J(0.4 V) of the p-Si/n-Si(χ) diode falls monotonically from 2.74909e-4 to 2.64981e-4 A/cm² for χ
    from 3.85 to 4.25 eV (legacy 2.749e-4 and 2.650e-4). The isotype barrier gives 6419.3, 4344.5 and 3114.4 A/cm²
    at 0 and ∓0.2 eV (legacy 6.42e3, 4.34e3, 3.11e3). Without band-gap narrowing the band shift is algebraically
    the legacy affinity gauge, so these are reproductions, not new numbers.
  - **Gate on a heterostructure:** the gate sits on silicon while node 0 is silicon with χ + 0.3 eV. At the
    silicon's flat-band voltage the gate charge is 3e-23 C/cm², where omitting the shift gives C_ox × 0.3 V =
    2.1e-7.
  - **Fermi–Dirac and band-gap narrowing:** n+ GaAs (1e19) / p-Al0.3Ga0.7As with thermionic emission. The electron
    Fermi level is flat to 8.9e-16 V at equilibrium, and the forward sweep to 0.6 V converges.
  - **Dimension:** a y-uniform 2D GaAs/AlGaAs diode reproduces 1D to 1e-9 at 1.2 V with thermionic emission.
  - **Bit-identity:** the Unit 14 hash probe (four silicon configurations) gives `main`'s hashes with the Newton
    floor (6.2, Unit 15) disabled. With it on, two hashes differ only in the recorded update values at −2 V. Every
    field, current and iteration count is identical to 17 digits.
- OLD / NEW / REASON:
  - **The band gauge.**
    - OLD: the default gauge was "nie", in which χ never reached the equations (a 0.5 eV step moved the solution by
      exactly 0). A "band_offset='affinity'" option gave the physical band edges with
      s = ln(Nc/n_ie,eff) + χ/V_T.
    - NEW: one gauge, always physical, with s taken from the un-narrowed depth.
    - REASON: the χ-blind gauge was kept only for legacy bit-identity. Built from the narrowed n_ie, the legacy
      affinity shift put all of the band-gap narrowing into the valence band and changed homojunction physics with
      the gauge. The un-narrowed depth keeps the narrowing split as in Unit 11, and it equals the legacy affinity
      gauge when narrowing is off.
  - **Compositions.**
    - OLD: the affinity gauge refused Fermi–Dirac (and incomplete ionization).
    - NEW: they compose.
    - REASON: the shift only moves η, so every statistics function applies unchanged.
  - **What counts as an interface.**
    - OLD: thermionic emission keyed interfaces on the material object, needed the affinity gauge, and refused a
      homojunction.
    - NEW: an interface is a parameter difference, and the model does nothing on a device of one material.
    - REASON: two equal parameter sets are one material; a model switch should not fail by device.
  - Kept from the legacy at first: field mobility on a heterointerface edge was refused. Superseded in the follow-up
    below.
- Not carried, or deferred:
  - Meshed insulator regions with a Si/SiO₂ interface, interface charge and traps, and interface recombination:
    track 15b, which the owner deferred (built in Unit 15b, 6.4).
  - The legacy's velocity-scaling knob, which its test used for the K → ∞ limit (a test hack, not an API).

**As built (Unit 15 follow-up, interfaces and outputs):**
- **Interfaces as device data (principle 9).** `device::Interface{region_a, region_b, transport}` with
  `InterfaceTransport::drift_diffusion` (the default for any edge) or `thermionic_emission`. `Device::create` rejects an
  unknown region, a region joined to itself, a pair declared twice (either order), a pair no mesh edge joins, and an
  unknown transport. `Device::transport(a, b)` looks a pair up.
- OLD / NEW / REASON, thermionic emission:
  - OLD: `PhysicsModels::thermionic_emission`, a global switch applying emission to every edge between two materials.
  - NEW: the transport is chosen per declared interface; undeclared material steps are drift-diffusion.
  - REASON: a graded composition described as many regions put emission on every step. With 10 steps the current fell to
    0.554 of the graded device's and depended on the step count. Now J(0.1 V) for 10, 20 and 40 steps is 33762.8,
    33901.3 and 33928.4 A/cm² (the 20-step value within 8e-4 of the 40-step one).
- Emission limits checked at device level: a declared interface without a step adds the series resistance V_T/(q v n)
  (measured 6.2575e-7 against 6.2654e-7 Ω cm²), and as A* grows (DOS, 1e4, 1e6, 1e8) J(TE)/J(DD) goes 0.034, 0.73,
  0.99947, 1.
- **Field mobility on a material step** (replacing the refusal): the edge mobility is the harmonic mean of each end's
  Canali mobility at the edge's field, each with its own material's parameters, with the exact Jacobian. On a
  one-material edge it is Unit 13's form.
- **Incomplete ionization** (`PhysicsModels::incomplete_ionization`, off by default as the legacy flag): the Poisson
  rows use N_D+ − N_A− at the electrons' and holes' own reduced energies (from n and p out of equilibrium, legacy), with
  the n, p (drift-diffusion) or ψ (equilibrium) derivatives; contacts and the neutral guess use
  `ionized_neutral_equilibrium`. The solved bulk's ionized fraction equals the 40-digit roots to 1e-10. A 4H-SiC p-n
  diode (aluminium 10% ionized) sweeps to 3 V.
- **Radiative recombination** (`PhysicsModels::radiative`, on by default since silicon's B is 0): R gains B(np − E) on
  physical densities, as Auger. A 200 µm-a-side GaAs diode with radiative recombination alone matches the long-base
  ideal-diode current with τ = 1/(B N) to 0.9% at 0.6–1.0 V, and with no recombination the short-base current q n_i²(D_n
  + D_p)/(N W) to 2.7e-4.
- **Band diagram** (`results::BandDiagram` in `EquilibriumResult` and `BiasPoint`; `assemble/band_edges.hpp`): E_c, E_v
  and the electron and hole quasi-Fermi levels in eV, measured from the equilibrium Fermi level, by the selected
  statistics; independent of which material node 0 is in. A contact at bias V holds its carriers' Fermi level at −V (to
  1e-9 eV).
- **Current resolution** (`BiasPoint::terminal_current_resolution`): a bound on how far the current through any cut can
  differ from the terminal current, the sum over the nodes of |F_n + F_p| plus 8 ε times the cancelling flux and
  recombination terms. In the Unit 11 1e19/1e18 silicon diode the measured edge-to-edge spread is 1.5e-8 A/cm² and the
  bound 8.9e-6, a factor 600. A first version summed only the contact's edges and fell below the spread (8.3e-9); this
  test caught it.
- `mesh::straddle_interface(axis, position, spacing)` puts two nodes at position ∓ spacing/2 and none between, since an
  interface lies at the midpoint of the edge joining two regions.
- Mutation checks, each caught (14): the radiative term or its dR/dn, either ionization Jacobian term, the degeneracy, a
  declared thermionic interface ignored, the Richardson constant ignored, the resolution without the residual, the hole
  Fermi-level sign, the mixed-edge Canali form, an interface without a joining edge, straddle_interface keeping nodes,
  the interface transport missing from the run record, the AlGaAs share ignored.

### 6.4 Contact description

A contact is a named set of boundary nodes with a type. First unit: ohmic, with the charge-neutral equilibrium boundary
value. Gate contacts on a lumped oxide came with Unit 12 (below); Schottky contacts are deferred. Contact boundary
conditions are residual and Jacobian rows, so contact handling lives in `assemble`, and terminal current is computed
from fluxes in `assemble`/`solve`. It is not an `analysis` quantity.

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
  anything), the warning for doping above 1e19 cm⁻³ under Boltzmann statistics (there is no warning channel; Fermi–Dirac
  statistics are available since Unit 14, `models.fermi_dirac`), and `graded_mesh` (not requested).

**As built (Unit 12, gate contacts and the MOS capacitor; owner request after Unit 11):**
- The legacy never meshes the oxide. Every MOS structure in it (`moscap.MOSCapacitor`, `Device2D/3D.add_gate`,
  `unstructured_dd3d`) treats the oxide as a lumped capacitor: a Robin condition on ψ at the silicon surface nodes. NiTCAD
  does the same. Meshed insulator regions are not built: they would need per-cell regions and permittivities in the mesh
  (the mesh has node regions only), and no legacy gate exercises them. They stay deferred.
- `ContactKind::gate` with a `GateStack`: the boundary patch the gate lies on (its face areas are the gate areas), t_ox,
  the oxide permittivity (3.9, the legacy `EPS_OX_R`), the electrode (n+ poly: φ_m = χ; p+ poly: χ + Eg(T); metal: a given
  φ_m, as the legacy `gate` argument) and the fixed oxide charge Q_f. Validation: the stack's numbers, the patch exists, and
  every node is on it. `SemiconductorParameters` gains `electron_affinity_eV` (silicon 4.05, legacy `chi`).
- The topology check now needs an *ohmic* contact in every connected part: a gate fixes no carrier density.
- `assemble/gate.hpp`: a gate node's Poisson row gains G_i(ψ_G,i − ψ_i) + S_i, with G_i the legacy κ times the node's
  scaled face area, ψ_G = (V_G − Φ)/V_T, Φ = φ_m − χ − (E_c − E_i) with E_c − E_i = Eg/2 + (kT/2) ln(Nc/Nv), and S_i the
  fixed charge on the face. The continuity rows keep zero boundary flux, and a gate carries no current. The gate charge
  (the oxide flux) is reported per bias point (`BiasPoint::gate_charge`, C/cm^(3−D)) and at equilibrium
  (`EquilibriumResult::gate_charge`).
  - Both assemblers share `GateNodes` (gate.hpp: per node the gate, its term and ψ_G; the row term; the charges).
  - Both assemblers' `set_bias` and the sweep's up-front check share `check_contact_bias` (contact_bias.hpp). The sweep's
    error names the point (context index), the contact (message) and the bias (value).
- OLD / NEW / REASON, intrinsic level (after the PR #16 review):
  - OLD: the legacy V_FB uses φ_semi = χ + Eg/2 − ψ_b V_T, i.e. it puts the intrinsic level at midgap.
  - NEW: ψ is referenced to the intrinsic level (n = n_ie e^ψ), which lies (kT/2) ln(Nv/Nc) above midgap. Φ includes
    that term, so every NiTCAD C-V and Id-Vg curve is the legacy one moved to higher gate voltage by 1.04 mV for
    silicon at 300 K. The legacy reproduction gate compares at V_G − 1.04 mV and still agrees to 1e-15.
  - REASON: the midgap form is inconsistent with the ψ both codes use. The error grows with the Nc/Nv asymmetry, e.g.
    about 39 mV for GaAs. Slotboom narrowing moves both band edges by half the narrowing, so it leaves E_i in place.
- OLD / NEW / REASON, gate reference:
  - OLD: the legacy row is κ(V_G − V_FB − (ψ_0 − ψ_b)). V_FB is computed for a substrate doping passed in separately,
    and ψ_b is the neutral potential of the node's own doping (`Device2D`, `psi_b_local`).
  - NEW: κ(V_G − Φ − ψ_0) plus the fixed charge. ψ_b cancels. Apart from the intrinsic-level term, the two forms are
    identical wherever the doping under the gate equals the substrate doping (every legacy MOS-C gate). The hand-written
    rows, with V_FB raised by the 1.04 mV, agree to 4.4e-16.
  - REASON: the electrode's Fermi level does not depend on the doping under it. The legacy form shifted the gate potential
    wherever the source and drain tails reach under the gate.
- OLD / NEW / REASON, gate charge:
  - OLD: the legacy `Qg` = C_ox(V_G − V_FB − φ_s) includes +qQ_f.
  - NEW: `gate_charge` is the charge on the electrode itself (the oxide flux), so Q_G + qQ_f + Q_semiconductor = 0.
  - REASON: it is the physical terminal charge. The capacitance dQ/dV is the same.
- `solve::Equations::equilibrium_poisson`, a `BiasOptions` choice: each point is thermal equilibrium (Poisson alone, every
  ohmic contact at 0 V, gates biased), as the legacy MOS-C solves it. That is the quasi-static C-V. `EquilibriumPoisson`
  gains `set_bias` (ohmic contacts must be at 0 V) and `gate_charges`; the run record lists `equations`. The sweep reads
  only the potential of an initial state, so it accepts one without densities.
- The run identity digests only what the sweep reads. In the quasi-static sweep it leaves out the mobility, SRH and Auger
  switches (and their settings) and the initial densities, and it never digests the work-function field of a
  polysilicon gate. Material parameters are digested whole.
- **Finding, a drift-diffusion MOS-C does not converge above threshold.** On the legacy MOS-C, Newton on (ψ, n, p) converges
  through accumulation and depletion, then stalls at V_G = 0 V, just below threshold (V_th = 0.093 V). The largest
  corrections sit on bulk minority electrons near the depletion edge (n ≈ 6e-14 scaled), and the residual of the surface
  electron row stays at 1e-4 to 1e-9.
  - Cause: the inversion layer has no ohmic contact of its own. Its electrons reach the substrate only through the depleted
    region, where n is some 1e-14 of the inversion density, so the electron subsystem is beyond double precision.
  - This is the physics of the low- versus high-frequency C-V, met as conditioning. The legacy never ran a MOS-C through
    drift-diffusion; its `Device2D` floors the relative density update at 1e-10 (M11-S5), which would not reach this
    (the corrections are 1e-3 relative on 1e-14 densities).
  - The quasi-static solve is the supported path for a MOS-C. A MOSFET, whose channel connects to the source and drain,
    converges through drift-diffusion in 5–7 iterations per point. A test records the stall so that a change is noticed.
- **Finding, subthreshold swing of the legacy MOSFET:** the legacy test records about 59.6 mV/decade for its MOSFET (5 nm
  oxide, N_A = 1e17). That is the thermal limit, with no body effect. NiTCAD gives 68.8 mV/decade on the same fixture,
  which is V_T ln 10 (1 + C_dep/C_ox) with C_dep at φ_s ≈ 0.7 V (68.1–71.7 between 2φ_F and φ_F). The legacy value is not
  consistent with its own parameters; the cause in the legacy code was not investigated (it cannot be run here). The
  legacy band (55–120) is kept and the body-factor bracket added.
- The MOSFET gates (about 9500 nodes, about 50 s in Release with Eigen SparseLU, 0.27 s per factorization) carry the hidden
  Catch2 tag `[.mosfet]` and run as the `solve_mosfet` ctest entry, in the Release configuration only.

**As built (Unit 15b, meshed insulators, electrodes and semiconductor-insulator interfaces; owner request, scope approved
by the owner with constraints on units, signs and exclusions):**
- `Region::material` is `std::variant<physics::Semiconductor, physics::Insulator>` (`device::Material`). `Device`
  gains `is_insulator(node)`, `relative_permittivity(node)`, `reference_node()` and `find_interface(a, b)`;
  `material(node)` requires a semiconductor node.
- `ContactKind::electrode` with `Electrode{kind, work_function_eV}`: electrostatic only (ψ fixed, no carrier rows, no
  current). Its nodes must be insulator nodes on a boundary patch. A polysilicon electrode is silicon's (φ_m = χ or
  χ + Eg of `silicon_parameters`): it has no semiconductor under it to take them from, unlike a lumped gate.
- `Interface` gains `fixed_charge_cm2` (Q_f, a sheet number density; the charge is q Q_f), `traps` and
  `recombination_velocity_n_cm_s` / `_p_cm_s`. They act on the semiconductor side of a semiconductor-insulator
  interface.
- Validation (each with a test): no semiconductor region; doping on an insulator node; an ohmic or gate contact node in
  an insulator, an electrode node in a semiconductor; an unknown electrode kind or a metal work function not positive;
  thermionic emission with an insulator side; charge, traps or recombination on an interface that is not
  semiconductor-insulator; Q_f not finite; a recombination velocity negative or not finite; traps outside the gap, with
  a density, cross-section or thermal velocity out of range.
- OLD / NEW / REASON, topology:
  - OLD: every connected part of the mesh graph needs an ohmic contact.
  - NEW: every connected part of the semiconductor nodes (joined by semiconductor-semiconductor edges) needs an ohmic
    contact, and every connected part of the graph an ohmic contact or an electrode.
  - REASON: carriers do not cross an insulator, so a silicon island enclosed by oxide has undetermined densities even
    when the oxide joins it to a contacted part; an oxide part needs its potential fixed. Without insulators the two
    rules are the old one.
- `results`: `interface_trap_charge` per declared interface in `EquilibriumResult` and `BiasPoint` (C/cm^(3−D), the
  traps' charge without Q_f); `gate_charge` includes electrodes; densities 0 and bands NaN on insulator nodes. The
  sweep's initial state need not have positive densities on insulator nodes. The run record hashes insulators,
  electrodes and interface charge, traps and recombination, only when present, so existing digests are unchanged.
- Bit identity: seven runs without insulators (1D and 2D silicon diodes, field mobility with Fermi–Dirac, a
  GaAs/AlGaAs thermionic junction with ionization, lumped MOS-C quasi-static in 1D and 2D and drift-diffusion) hash every
  field, band, current, charge, convergence record and run digest identically on `main` and on the branch.
- **Finding, legacy D_it charge:** `moscap.cv_sweep` with D_it reports Q_g = C_ox (V_G − V_FB − φ_s) − q D_it φ_s. By
  Gauss, C_ox (V_G − V_FB − φ_s) is already the gate charge and the trapped charge is −q D_it φ_s, so the legacy Q_g is
  −Q_s, the semiconductor charge, and its C = dQ_g/dV_G leaves out the traps' response, which a measured quasi-static
  C-V contains. NiTCAD reports the electrode's own charge. The legacy's φ_s(V_G), which the trap term does shape
  correctly, is the one reproduced (6.2, Unit 15b).

### 6.5 Input representation

The authoritative input is the in-memory C++ `device` description plus a solver-options struct. An
on-disk input format is **deferred** (R2) until concrete requirements exist. The input types must not depend on
any file format: no serialization types or format concepts inside them, and any future reader or writer depends on
them, never the reverse. The legacy wire format was a Python JSON `DeviceSpec` and is not carried over.

### 6.6 Result representation

`results` holds fields indexed by mesh node, terminal quantities per bias point, and a run record
(inputs hash/identity, options, per-iteration convergence history, converged flag). It is plain data.
An on-disk result format is **deferred** (R2), under the same independence rule as 6.5.

**As built (Unit 10, `include/NiTCAD/results/`, header-only, depends on `base` only):**
- `NodeFields`: potential, n and p per node.
- `IterationRecord` and `ConvergenceRecord`: per iteration the full correction and the residual before it; the converged
  flag; the smallest pivot ratio.
- `EquilibriumResult`: fields and convergence.
- `BiasPoint`: the bias per contact, fields, terminal currents, edge currents and convergence.
- `RunRecord`:
  - `input_identity`, a 64-bit FNV-1a digest of the mesh, device, options, bias points and initial state, hashing
    doubles by their bit pattern. It identifies a run on one build; it is not a cross-machine checksum (6.8).
  - `settings`, the options as (name, value) pairs. These are plain pairs because `results` may not depend on `solve`.
- `Sweep`: the run record, the completed points, and `stopped` (the error that ended the run early, if any) with
  `unfinished` (the convergence history of the point being solved then).
- No file format (R2). `results` has no test directory of its own: it is plain data, exercised through the solve tests.
- Later units add fields: the band diagram and current resolution (Unit 15), the interface trap charge (Unit 15b).
- `Transient` (Unit 21, `results/transient.hpp`): the run record, a `TimePoint` per accepted step (time, step, order,
  error ratio, bias, total, conduction and displacement current per contact, contact charges, trapped charges,
  convergence), `TransientSnapshot`s (fields, band diagram, edge currents, trap occupancies) at t = 0, the output times
  and t_end (or every step on request; every step of a large run would be gigabytes), the rejected-step count, and
  `stopped` / `unfinished` as `Sweep`.

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

**As built (Unit 10, `solve/control.hpp`, `solve/bias.hpp`):**
- `RunControl{std::stop_token stop; std::function<void(const Progress&)> progress}`. It is accepted by `solve_equilibrium`,
  `solve_bias` and `sweep_bias`.
- **Safe points:** the token is checked before every Newton iteration, the first included, and between bias points.
- **Progress:** an event after every Newton iteration, with phase (equilibrium or bias), point, point count, iteration,
  update, residual and converged flag. Events are strictly increasing in (phase, point, iteration), and the bias events are
  exactly the points' convergence records. The callback runs on the solving thread; it must not throw or re-enter the solve.
- **`sweep_bias(device, points, options, initial, control)`:**
  - It solves the points in order, each started from the previous one, with one system and one linear solver, so the
    pattern is analyzed once.
  - Invalid input is an `std::expected` error before anything runs.
  - Once running, a cancelled or failed sweep is still a value: it keeps the completed points, with `stopped` (`cancelled`
    and the point or iteration in the context, or `non_convergence`, `singular_system`, ...) and `unfinished`.
  - `solve_bias` is a sweep of one point and returns a stopped run as its error. `newton_solve` fills a
    `ConvergenceRecord` whatever the outcome, and takes the stop token and an iteration observer.
- Not built here: the worker thread and the Win32 message posting belong to `app` (L8). A test runs a sweep on a
  `std::jthread` and cancels it from the main thread deterministically (the worker pauses in a progress event until the
  stop is requested).
- Cancellation latency is at most one residual evaluation, factorization and solve, as stated above; the factorization
  times measured in 6.10 (Known limits) bound it.
- Unit 21: `solve_transient` takes the same `RunControl`. The token is also checked before every step attempt; a
  cancelled run keeps its accepted steps. `Phase::transient` events carry the step attempt as point (point count 0,
  not known in advance) and the step's end time in the new `Progress::time_s`.

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

### 6.12 Long-term computational architecture

The following are **target capabilities**, not currently implemented layers. They are recorded here so future units
do not accidentally make the core device-specific.

#### 6.12.1 Generic fields

A future field abstraction should describe at least: name, units, location (node/edge/face/cell), scalar/vector/tensor
shape, known versus unknown status, region applicability, and storage. Specialized `NodeFields` and solution types
remain valid until a generic representation is actually needed.

#### 6.12.2 Generic equations

A future equation abstraction should identify its unknown fields, domain, residual contribution, Jacobian contribution,
and any boundary/interface terms. The assembler may then activate a set of equations from the device/model description
and build the same global `F(x)` and `J(x)` machinery used by the current drift-diffusion system.

The intended consequence is that these are compositions rather than separate solvers:

```text
Diode     = Poisson + electron continuity + hole continuity
MOS       = diode equations + oxide Poisson + gate/interface conditions
MOSFET    = MOS equations + multiple contacts + transport + requested mobility/statistics models
Thermal   = existing semiconductor equations + heat equation + thermal material/BC models
Hydrodynamic = existing equations + carrier-energy equations + corresponding constitutive models
```

The exact equation sets remain deferred until their units are requested.

#### 6.12.3 Automatic or assisted differentiation

The current contract remains explicit value + exact partial derivatives. Long-term, NiTCAD may support automatic or
symbolic differentiation for complex model composition while retaining hand-written kernels on measured hot paths.
Any such mechanism must preserve the existing Jacobian accuracy contract and must not leak into the public layer
dependencies unless a unit requires it.

#### 6.12.4 Adaptive mesh refinement

The existing graph-only mesh interface is the foundation for future unstructured and adaptive meshes. The target loop
is:

```text
solve → estimate error / feature indicator → mark → refine/coarsen → transfer state → solve
```

Refinement should be able to concentrate resolution around junctions, interfaces, high-field regions, tunnelling
regions, impact-ionization regions and other requested features without changing the physics or solver interfaces. Mesh
transfer must preserve the physical units and the conservation/continuity gates appropriate to the active equations.

#### 6.12.5 Solver hierarchy

The current Newton and backend-neutral linear-solver contracts are the base. Future solver strategies may include damped
Newton, line search, Gummel, pseudo-transient continuation, homotopy/continuation, iterative linear solvers, algebraic
multigrid and block preconditioners. They are alternatives behind solver interfaces, not device-specific algorithms.

#### 6.12.6 Verification and conservation

A future verification layer should report more than `converged`: nonlinear residual, linear-solve accuracy, Jacobian
consistency, current/charge conservation, relevant energy balance, mesh quality and other physics-specific invariants.
It belongs above numerics conceptually and should consume results rather than modifying the solver.

#### 6.12.7 Scientific results and provenance

The current plain-data `results` design should grow toward a scientific data model containing fields, terminal
quantities, sweeps, time/frequency axes when applicable, solver diagnostics, convergence history and provenance. On-disk
formats remain decoupled from internal types (R2).

#### 6.12.8 Device neutrality rule

No future unit may introduce classes whose primary purpose is to hard-code a device family (for example `MOSFETSolver`
or `HEMTSolver`) when the capability can be represented by the existing region/interface/field/equation/model concepts.
A device family may have a convenience builder, example or analysis workflow, but the numerical engine must remain
composable.

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
include/NiTCAD/physics/     semiconductor, mobility, recombination, statistics (Unit 5); bandgap_narrowing (Unit 11);
                            electron affinity (Unit 12); field_mobility (Unit 13); fermi_dirac, Fermi-Dirac
                            statistics (Unit 14); material sets, thermionic_emission, ionization (Unit 15)
src/physics/                parameter validation, band quantities, models; increasing_root (private) (Units 5, 11,
                            13, 14, 15)
tests/physics/              published values, limits, FD derivative gates; heavy doping; Canali; Fermi integral and
                            Fermi-Dirac statistics; materials and emission velocity (Units 5, 11, 13, 14, 15)
include/NiTCAD/device/      device.hpp, contact.hpp (Unit 6); gate contacts (Unit 12)
src/device/                 description validation, topology check           (Unit 6, exists)
tests/device/               construction, validation, floating regions       (Unit 6, exists)
include/NiTCAD/assemble/    scaling, bernoulli, sg_flux, ohmic, equilibrium_poisson (Unit 7); drift_diffusion (Unit 9);
                            models (Unit 11); gate, contact_bias (Unit 12); thermionic_flux, band_edges (Unit 15)
src/assemble/               scaling, scaled device (private), equilibrium Poisson, drift-diffusion, gate, contact bias
                            (Units 7, 9, 12)
tests/assemble/             Bernoulli/SG references, FD-Jacobian gate, reduction (Unit 7); gate terms (Unit 12);
                            heterojunctions (Unit 15)
include/NiTCAD/solve/       newton.hpp, equilibrium.hpp (Unit 8); bias.hpp (Unit 9, sweeps Unit 10); control.hpp (Unit 10)
src/solve/                  equilibrium and bias solves, sweeps, run record (Units 8-10); fields helper (private,
                            Unit 12)
tests/solve/                Newton contract, equilibrium and bias diode gates, legacy graded_mesh port, sweeps,
                            cancellation and progress (Units 8-10); MOS-C (legacy moscap port) and MOSFET (Unit 12);
                            field mobility (Unit 13); Fermi-Dirac (Unit 14); heterojunctions (Unit 15)
include/NiTCAD/results/     convergence.hpp, solution.hpp, run.hpp (header-only plain data) (Unit 10, exists)
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

## 11. Build order and long-term capability roadmap

Each implemented unit is one branch, is independently tested, and is merged only on the owner's instruction. The table
below records the completed foundation and the **target capability order** for future units. Future ordering is a design
direction, not an automatic permission to implement every row.

Unit 1 was named by the owner on 2026-10-02 and has no open prerequisites (V1–V4 are done; V3 is consumed at Unit 7). The
order is by dependency, so each unit is
testable on arrival. Units 4–9 are dimension-generic from the start (D4) and use no tensor-grid assumptions (6.11).

| # | Unit | Layer | Legacy reference | Gate |
|---|---|---|---|---|
| 1 | Build scaffold: CMake, C++23 gate, vcpkg manifest with pinned baseline, Catch2 v3 harness (**done, on `main`**) | — | `core/CMakeLists.txt` (reference for options only) | configure fails on a toolset without the needed C++23 features; dependencies resolve from the pinned manifest; one trivial test runs |
| 2 | Constants, units, `Error` type, `NITCAD_EXPECTS` (**done, on `main`**) | base | `constants.py`, `core/include/tcad/base/errors.hpp` | CODATA 2018 values; V_T(300 K) = 0.025852 V |
| 3 | Sparse matrix, backend-neutral solver interface (6.10), Eigen SparseLU backend (**done, on `main`, with the hardening**) | linalg | `linsolve.py`, `core/src/solver/direct_lu.cpp` | known systems (e.g. analytic tridiagonal); accuracy check (backward error after hardening); singular system returns an error; symbolic-reuse path exercised; header-boundary check, no second backend (Q5) |
| 4 | Generic node/edge/control-volume mesh; tensor-grid constructors for D = 1, 2, 3 (**done, on `main`**) | mesh | `mesh.py`, `mesh2d.py`, `mesh3d.py`, `core/include/tcad/mesh/stencil.hpp` | total volume matches domain (section 10); positive dual volumes; consistent edge geometry across D = 1, 2, 3 on a uniform grid; public interface contains no (i, j, k) indexing |
| 5 | Si material parameters, Caughey–Thomas mobility, SRH recombination, Boltzmann statistics (value + partials) (**done, on `main`**) | physics | `materials.py`, `core/include/tcad/physics/` | n_i(300 K) ≈ 1.0674e10; published mobility values (legacy `test_caughey_thomas_matches_published_silicon_values`); derivatives vs finite differences; SRH vanishes at equilibrium |
| 6 | Device description and ohmic contact data (**done, on `main`**) | device | `core/include/tcad/device1d/inputs.hpp` (reference only) | construction and validation; invalid input returns an error |
| 7 | Scaling (V3), SG/Bernoulli flux, equilibrium Poisson residual + Jacobian, ohmic boundary (**done, on `main`**) | assemble | `core/src/device1d/device1d.cpp`, `core/src/device1d/inputs.cpp` | Bernoulli limits and symmetry; FD-Jacobian gate (section 10) |
| 8 | Newton solver (scaled variables) and equilibrium solve (**done, on `main`**) | solve | `device.py` options, `device1d.cpp` | built-in potential within 2e-3 V; bulk neutrality; convergence; non-convergence returns an error value |
| 9 | Electron/hole continuity assembly and bias solve; first end-to-end gate (**done, on `main`**) | assemble, solve | `device1d.cpp`, `tests/test_validation.py`, `tests/test_device1d_native_gates.py` | J(0.5 V) = 1.280e-2 A/cm² ± 1%; ideal-diode law; current continuity; mesh independence; uniform 2D/3D reproduces 1D to a tolerance set at this unit |
| 10 | Result representation, cancellation and progress (**done, on `main`**) | results, solve | — (new) | cancellation returns partial results; progress is monotonic |
| 11 | Slotboom band-gap narrowing and Auger recombination (**done, on `main`**) | physics, assemble | `materials.py` (`bandgap_narrowing_slotboom`, `nie_effective`, `recombination`), `device1d.cpp` (`delta`, `delta_p`) | legacy model benchmarks; Slotboom and n_ie against 40-digit values; Auger partials vs FD; zero equilibrium current with varying n_ie; diffusion current × exp(ΔEg/kT); FD-Jacobian gate on a 1e19/1e18 diode |
| 12 | Lumped-oxide MOS: gate contacts on a lumped oxide (no meshed oxide region), quasi-static (equilibrium) sweeps, the MOS capacitor and MOSFET (**done, on `main`**, `1de5967`) | physics, device, assemble, solve | `moscap.py` (`MOSCapacitor`, `flatband_voltage`), `device2d.py` (`GateBC`, `add_gate`), `mosfet.py` | legacy C-V physics checks P1–P9, temperature and fixed-charge checks; NiTCAD reproduces the ported legacy MOS-C solve; legacy MOSFET threshold, on/off, swing, monotonic and mesh-independence gates |

**Target capabilities (Units 13–27).** Design direction, not permission to build: each needs an owner-named
unit, and the order may change (11.1). Meshed oxide regions and a Si/SiO₂ interface, beyond Unit 12's lumped
oxide, belong to the heterojunction and interface track (15); the owner split them into a later unit (15b), now
built.

| # | Unit / capability | Layer(s) | Status / purpose |
|---|---|---|---|
| 13 | Field-dependent mobility / velocity saturation | physics, assemble | **done, on `main`** (`f8bda70`): Canali per edge with exact Jacobian (6.2, Unit 13) |
| 14 | Fermi–Dirac statistics and high-density carrier models | physics, assemble | **done, on `main`** (`8610d9a`): parabolic-band Fermi–Dirac statistics, the legacy ν-factor scheme with an exact Jacobian (5 and 6.2, Unit 14); incomplete ionization deferred (14.3) |
| 15 | Heterojunctions, band offsets and interface transport | device, assemble, physics | **done, on `main`** (`2a8718a`): band offsets through a per-node band shift, permittivity steps, thermionic emission, the legacy material sets (6.3 and 6.2, Unit 15); follow-up: interfaces as device data, incomplete ionization, radiative recombination, band diagram, current resolution; meshed insulators, interface charge and traps are 15b |
| 15b | Meshed insulators, semiconductor-insulator interfaces, interface charge, traps and recombination | physics, device, assemble, solve | **done, on `main`** (`4e23119`): insulator regions, electrodes on a meshed oxide, fixed charge, interface traps (levels and uniform bands, steady-state SRH occupancy), surface recombination, all at the interface potential (5, 6.2 and 6.4, Unit 15b) |
| 16 | Unstructured mesh | mesh, assemble | target |
| 17 | Adaptive mesh refinement and state transfer | mesh, solve, results | target |
| 18 | Scalable linear-solver backends: PARDISO and/or iterative/AMG paths | linalg | target; must preserve backend-neutral interface |
| 19 | Impact ionization and breakdown-oriented continuation | physics, assemble, solve | target |
| 20 | Band-to-band tunnelling and nonlocal path machinery | assemble, physics, solve | target |
| 21 | Transient simulation | solve, assemble, results | **done on branch `solve/transient`**: backward Euler and variable-step BDF2 with error-controlled steps, waveforms, displacement current with exact conservation, interface trap dynamics eliminated in the interface solve (6.2, Unit 21) |
| 22 | AC small-signal analysis | linalg, assemble, solve, results | target; complex system already reserved by A5 |
| 23 | Thermal / electrothermal coupling | physics, assemble, solve | target |
| 24 | Analysis and extraction engine | analysis | target |
| 25 | Scientific visualization / rendering | render | target |
| 26 | Native Windows application and workflow | app | target; Win32 API, x64 target, Direct3D 12 / Direct2D |
| 27 | Optimization / sensitivity / inverse-design workflows | analysis, solve | long-term target |

Units 1–9 are the smallest end-to-end vertical slice: a validated drift-diffusion diode. Everything past
Unit 10 is deferred until sequenced. Units 11–15, 15b and 21 were requested by the owner after Unit 10; nothing
else past them is started.

### 11.1 Capability tracks may interleave

The numerical architecture should not be forced into one linear feature list. A future unit may be pulled forward when
it unlocks a critical downstream capability, but it must preserve these priorities:

1. strengthen the generic equation/field/interface abstractions before multiplying device-specific special cases;
2. strengthen numerical robustness and performance before adding highly stiff physics;
3. keep mesh, physics, solver and application dependencies one-way;
4. require physics-specific verification rather than accepting convergence alone;
5. maintain dimensional generality and avoid a separate 1D/2D/3D implementation path.

### 11.2 North-star milestone

A mature NiTCAD simulation should look conceptually like:

```text
Device = regions + interfaces + contacts + materials + parameters
        ↓
Fields = unknown + fixed + derived quantities
        ↓
Equations = selected physics + boundary/interface laws
        ↓
Mesh = structured or unstructured, optionally adaptive
        ↓
Assembly = residual + exact/automatic Jacobian
        ↓
Solve = nonlinear strategy + linear backend + continuation
        ↓
Verification = numerical + conservation + physical checks
        ↓
Results = fields + terminals + sweeps + diagnostics + provenance
        ↓
Analysis / visualization / optimization
```

No device family should require a second numerical architecture to participate in this flow.

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
withdrawn; numbering is otherwise kept stable. This document is maintained against the current `main`
architecture; historical branch names remain only where they are useful to explain how a decision was reached.

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
| R4 | Boltzmann carrier statistics initially; Fermi–Dirac for parabolic bands added at Unit 14 on the owner's request (section 5, item 9); other statistics deferred until requested. |
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
| — | Statistics beyond parabolic-band Fermi–Dirac (Unit 14); incomplete ionization | Requested (R4). |
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
| V17 | Unit 10 (`solve/results-control`): Debug and Release build with no warnings; `nitcad_solve_test` has 32 test cases, all pass (the existing Unit 8 and 9 tests moved to the `results` types unchanged in substance). A 0–0.5 V sweep equals the chain of one-point solves to 1e-9 (the chain passes states through physical units) and ends at J(0.5 V) within the 1% gate; the run record carries the settings. Progress: strictly increasing in (phase, point, iteration), iterations counting from 1 per point, the converged event closing each point, the bias events equal to the convergence records. Cancellation keeps partial results: a stop after point 1 returns points 0–1 bit-identical to the full run (context: point 2, no unfinished point); a stop at iteration 2 of point 3 returns points 0–2 and the unfinished 2-iteration history (context: iteration 3); a stop requested before the start solves nothing and emits no event, and `solve_equilibrium` and `solve_bias` return `cancelled`; a sweep on a `std::jthread` cancelled from the main thread stops at the next safe point, deterministically. A point that does not converge (0.8 V with 8 iterations) stops the sweep with the earlier point kept. Invalid points or initial state are errors before any event. The run identity is equal for equal inputs and changes with the doping, an option, a bias point or an initial state. Mutation checks: removing the between-points or the Newton stop check, a wrong point index in progress events, and leaving the donors out of the digest each fail the suite. | Verified locally. |
| V18 | Unit 11 (`physics/bgn-auger`): Debug and Release build with no warnings; `nitcad_physics_test` 38 test cases, `nitcad_assemble_test` 24, `nitcad_solve_test` 35, all pass. Slotboom ΔEg at 1e18/1e19/1e20 and the effective n_ie match 40-digit values to 1e-13 and 1e-12; zero at and below N0; the legacy benchmarks pass (BGN positive and monotonic; SRH + Auger more than doubles when the densities double). Auger: zero at equilibrium, exactly ×8 for doubled high-injection densities, partials vs FD to 1e-8 with and without a carrier-dependent equilibrium product; the Auger part of the assembled Jacobian on its own vs FD (gate 1e-5); the coupled FD-Jacobian gate on a 1e19/1e18 diode with all models. Device level: equilibrium edge currents of a 1e19/1e18 diode within 0.0095 of their rounding bound (without the ln n_ie term the junction edges would carry about 1e8 A/cm²), V_bi from the effective n_ie to 1e-12; the 1e19/1e19 diffusion current ratio 8.5908 against exp(ΔEg/kT) = 8.5932 (gate 0.5%); Auger on the legacy fixture +4.3e-7 with J(0.5 V) still within 1%. The legacy rows written out by hand now include BGN (Δln n_ie) and Auger, on a 1e18/1e16 diode, to 1e-12. Mutation checks: no Δln n_ie term, Auger dR/dn dropped from the Jacobian, and the BGN flag ignored each fail the suite. | Verified locally. |
| V19 | Unit 12 (`device/mos`): Debug and Release build with no warnings; `nitcad_physics_test` 38 test cases, `nitcad_device_test` 18, `nitcad_assemble_test` 30, `nitcad_solve_test` 52 (+2 `[.mosfet]` cases, Release only), all pass. Gate terms: the scaled coupling equals the legacy κ times the scaled face area in 1D/2D/3D to 1e-15; the 1D rows of both assemblers equal the legacy MOS-C rows written out by hand (with V_FB and Q_f) to 4.4e-16; only the gate node's Poisson row changes and no carrier flux crosses the gate (bitwise); FD-Jacobian gate with gates, worst 6.7e-10 (1D), 2.2e-8 (2D), 3.3e-8 (3D), 2.1e-9 (equilibrium Poisson). MOS-C on the legacy fixture (1e17, 5 nm, n+ poly, 1200 nodes, −2 to 2 V): the quasi-static sweep reproduces a C++ port of the legacy solve to 5.6e-16 V in φ_s, 6.1e-16 V in ψ at every node, and 7.2e-15 C_ox in C; the legacy checks P1 (residual 4.9e-12, gate 1e-6), P2 (Gauss to 1.4e-14 relative, legacy 2%), P3 (1.1e-3 and 1.2e-2, gate 5%), P4, P5, P6 (C_min +10.1%, W_max −10.5%, gate 15%), P7 (flatband crossing 1.0 mV and 2φ_F crossing 0.35 mV from the landmarks, gate 50 mV), P8, P9 (a tighter Newton tolerance moves φ_s by 0) pass; C_min within 15% at 275, 300, 350 K; Q_f and the electrode work function shift the curve rigidly to 1e-16 V; y-uniform 2D and 3D reproduce 1D to 1e-12. Drift-diffusion equals the quasi-static state in accumulation and depletion (ψ to 1e-15 V, p to 2e-14) and stalls once inversion starts (recorded). MOSFET (legacy fixture, 9490 nodes, drift-diffusion): V_th by max g_m 0.151 V against the landmark 0.093 V (gate 0.1 V), on/off 2e11 (gate 1e6; the off current is at the 2e-12 A/cm rounding floor), swing 68.9 mV/decade (legacy band 55–120; body-factor bracket 68.1–71.7), monotonic from 0.3 to 1.5 V, Kirchhoff to 5e-12 A/cm, Id(1 V) on a 2× finer mesh within 0.055% (gate 10%). Mutation checks, each caught: the gate Jacobian term dropped (either assembler), the fixed-charge sign, an L_D^(D−2) area scale, p+ poly without Eg, the midgap offset without Eg/2, gates anchoring the topology, ohmic bias accepted at equilibrium, the gate bias ignored at equilibrium, the gate-charge sign, the patch check. | Verified locally. |
| V19a | Unit 12 after the PR #16 review (all ten findings applied): `nitcad_assemble_test` 31 test cases, `nitcad_solve_test` 55 (+2 `[.mosfet]`), Debug and Release with no warnings, all pass. Φ is now referenced to the intrinsic level: the offset matches φ_m − χ − Eg/2 + (kT/2) ln(Nv/Nc) to 1e-15 at 300 and 400 K, with a shift of 1.0416 mV for silicon at 300 K. The legacy reproduction holds at V_G − 1.0416 mV: φ_s to 6.7e-16 V, ψ to 4.2e-16 V, C to 5.1e-15 C_ox. P1 gives 6.0e-12 at the shifted biases; the flatband and 2φ_F crossings move by +1.04 mV (measured 0.05 and 1.34 mV from the landmarks; gate 50 mV), C_min +10.2%; the MOSFET V_th is 0.152 V and the swing 68.8 mV/decade. The drift-diffusion gate diagonal is accumulated (`-=`). Both assemblers use `GateNodes` and `check_contact_bias`; the edges and the gate share `permittivity_ratio`; `solve_equilibrium` and the quasi-static sweep share one fields helper. New tests: the contact-bias rule; the sweep error names point, contact and bias; a quasi-static warm start from a potential alone (drift-diffusion still rejects it); the run identity ignores the transport switches and initial densities in the quasi-static sweep and a polysilicon gate's work-function field, but not BGN or a metal work function; `solve_equilibrium` reports the zero-bias gate charge, equal to a one-point quasi-static sweep to 1e-12. Mutation checks, each caught: the intrinsic-level term dropped, the gate diagonal dropped, ohmic bias accepted at equilibrium, densities required for a quasi-static start, the SRH switch digested in the quasi-static sweep, the poly work function digested, the equilibrium gate charge missing, the contact missing from the sweep error. | Verified locally. |
| V20 | Unit 13 (`physics/field-mobility`): Debug and Release build with no warnings; `nitcad_physics_test` 44 test cases, `nitcad_assemble_test` 34, `nitcad_solve_test` 61 (+3 `[.mosfet]`, Release only), all pass. Canali: the legacy formula for β = 1, 2 and 1.5 to 1e-15; μ(0) = μ0, μ falling and μE rising monotonically to v_sat (within 1e-3 at 1e9 V/cm); the low-field expansion; dμ/dE against fourth-order differences to 4.2e-12 of μ/E; the right-hand derivative at E = 0; v_sat and β validated. Assembly: FD-Jacobian gate with field mobility 2.7e-9 (1D), 2.3e-9 (2D), 6.0e-9 (3D) (gate 5e-5); the field-mobility part on its own (ψ columns) 2.5e-7 (gate 1e-5); every edge's current scaled by exactly its Canali factor (2.2e-16; factors down to 0.063 at the probe). Device level: a uniform n- and p-type resistor (100 nm, 1e16, 0.1 V steps to 1e6 V/cm) matches q(nμ_n(E) + pμ_p(E))E to 5.8e-15; electrons reach 0.99996 v_sat and holes 0.9798 v_sat (x/(1 + x) for β = 1); y-uniform 2D and 3D reproduce 1D to 1e-12; at most 5 (electrons) or 7 (holes) Newton iterations per point against 2 with constant mobility; quadratic finish (2.8e-2 then 1.2e-8). Legacy diode: J(0.5 V) −0.77% with the model on, total current the same on every edge to 7.5e-8. MOSFET at V_G = 1 V (max_update 1): I_D ratio 0.968 at 0.1 V falling monotonically to 0.735 at 1 V. Mutation checks, each caught: the field-derivative term dropped, a wrong field scale, the Canali derivative without x^(β−1), the switch ignored, β < 1 accepted, the switch missing from the run record, the hole parameters taken from the electrons. | Verified locally. |
| V21 | Unit 14 (`physics/fermi-dirac`): Debug and Release build with no warnings; `nitcad_physics_test` 56 test cases, `nitcad_assemble_test` 40, `nitcad_solve_test` 69 (+3 `[.mosfet]`, Release only), all pass. F₁/₂ within 1.5e-15 of 40-digit references from −700 to 1000; F₋₁/₂ within 1.3e-13 to η = 15 and 3e-12 to 40; the returned derivative matches central differences to 1e-8 across both joins; ln F₁/₂ increasing and concave to rounding on a 1e-3 grid over −50 to 100; the inverse to 1e-14; Boltzmann and Sommerfeld limits with their published correction terms. Statistics: densities, degeneracy factors and the equilibrium product against finite differences; generalized mass action to 1e-13; the neutral equilibrium against 40-digit roots to 1e-13 (η) and 1e-12 (n, p). Device gates and mutation checks in 6.2 "As built (Unit 14)". The Boltzmann path is bit-identical to `main` (probe hashes, 6.2). | Verified locally. |
| V22 | Unit 15 (`device/heterojunctions`): Debug and Release build with no warnings; `nitcad_physics_test` 60 test cases, `nitcad_assemble_test` 47, `nitcad_solve_test` 77 (+3 `[.mosfet]`, Release only), all pass. Material sets, depths and the emission velocity against 40-digit values; the abrupt-heterojunction first integral (V_bi to 1e-12, interface potential within 1.4e-5 V, D within 7.6e-4); legacy M33 G1, G2 and S2 gates, with G2 currents reproduced to 4 digits; FD-Jacobian gates with interfaces, thermionic emission and Fermi–Dirac in 1D/2D/3D; the gate flat band on a heterostructure; 2D = 1D. Details, the Newton floor and the ten mutation checks in 6.2 and 6.3 "As built (Unit 15)". | Verified locally. |
| V22a | Unit 15 follow-up (`device/heterojunctions`): Debug and Release build with no warnings, all suites pass (Release 9/9 with `solve_mosfet`, Debug 8/8). FD-Jacobian gates with incomplete ionization (Poisson 1.9e-9; drift-diffusion 1D/2D/3D at most 1e-8), the ionization part on its own 1.4e-7 and 2.1e-7 (gate 1e-5), the radiative part 4.8e-9, thermionic part 2.0e-9, a non-planar 2D interface. Ionized fractions against 40-digit roots to 1e-10; radiative long-base GaAs diode within 0.9% of the analytic current, short-base without recombination to 2.7e-4; emission resistance 6.2575e-7 vs 6.2654e-7 Ω cm²; J(TE)/J(DD) → 1 as A* grows; graded staircase 10/20/40 steps converging (8e-4); current resolution bounds the spread (factor 600); band diagram Fermi levels at the contacts' biases to 1e-9 eV. Fourteen mutation checks caught. | Verified locally. |
| V23 | Unit 15b (`device/insulators`): Debug and Release build with no warnings; `nitcad_physics_test` 72 test cases, `nitcad_device_test` 24, `nitcad_assemble_test` 57, `nitcad_solve_test` 96 (+4 `[.mosfet]`, Release only), all pass. Without insulators every output of seven probe runs hashes identically to `main`. Fermi occupancy, trap-band quadrature and Gauss–Legendre nodes against 40-digit values (4 ε; 7.1e-16; 25 digits); SRH occupancy equals the Fermi function at equilibrium within 1.9e-15; FD Jacobians 1D/2D/3D at most 7.0e-8, the interface part 2.4e-9; meshed MOS-C second order against the exact solution (7.6e-6 V at 2400 nodes); the Q_f shift equals its discrete value to 6e-14; the legacy D_it stretch-out within 0.28 mV; surface recombination against the analytic diode within 8.3e-4; meshed against lumped MOSFET thresholds 0.21 mV apart, currents within 0.58%. Seventeen mutation checks caught (the run-identity test first missed the electron thermal velocity and was extended). | Verified locally. |
| V23a | Unit 15b interface potential (`device/insulators`, on `daf4de7`): Debug and Release build with no warnings, all suites pass. Q_f, traps and surface recombination at ψ_I by a local solve per interface edge with the exact Jacobian: D_it MOS-C second order (2.94e-4 to 4.39e-6 V, ratio 4.04–4.08, was first order), Q_f flat-band shift = −q Q_f/C_ox to 6e-13, FD Jacobians at most 1.7e-7 (interface part 2.0e-8), surface recombination within 8.3e-4, unchanged devices bit-identical to `main` and `daf4de7`, solve time within 0–20% (noise 10–20%). Twenty-six mutation checks caught. | Verified locally. |
| V24 | Unit 21 (`solve/transient`): Debug and Release build with no warnings; `nitcad_physics_test` 72 test cases, `nitcad_device_test` 24, `nitcad_assemble_test` 63, `nitcad_solve_test` 112 (+5 `[.mosfet]` and 5 `[.transient]`, Release only), all pass (Release 10/10, Debug 8/8). Time-step FD Jacobians 1D/2D/3D, both statistics, incomplete ionization, at most 1.7e-7 (trap part 3.5e-7); a long step is the steady state; total currents sum to 3.1e-10 (2D MOS with traps and two ohmic contacts) and integrate to the contact charges; backward Euler first order (1.92–1.99 per halving) and BDF2 second (3.33 to 3.90); error control follows rtol across kinks and a jump; MOS-C RC response within 0.093% (meshed and lumped), dielectric relaxation within 0.34%, trap emission rate within 3.2e-4; the legacy diode turn-off (legacy core built from the reference checkout) within 3.7e-5 until the legacy run stalls (a legacy finding); 2D/3D extrusions to 2.3e-13; the meshed MOSFET gate step settles to the Id–Vg current within 1e-6. Steady paths bit-identical to `main` (plain, meshed-oxide and trap devices, run digests included). Twenty-seven mutation checks caught (four after a test was added or sharpened). | Verified locally. |
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

### 14.6 North-star architectural rule

The following is now the governing long-term design principle, subject to unit-by-unit verification:

> **Build NiTCAD as one composable semiconductor multiphysics equation engine. Do not build a collection of
> device-specific numerical solvers. Regions, interfaces, contacts, fields, equations, constitutive models, mesh,
> nonlinear solving, linear solving, verification and results remain reusable across device families and dimensions.**

This principle does not override the existing feature-freeze rule. An implementation change still requires an
owner-named unit and must be validated against the gates in this document.

**Approval workflow:** the owner names each unit. Each unit is built on its own branch, tested, committed and reported, and is
not merged without the owner's instruction.
