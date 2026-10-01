# NiTCAD — Architecture

**PROPOSAL, accepted as the current proposal. Decisions recorded (section 14). Verification V1–V3 was run on 2026-10-02 (section 14.4); its findings need owner decisions (Q6–Q9) before Unit 1.**

The owner accepted this architecture as the current proposal and recorded decisions D1–D6, D8 and R1–R5.
Section 14 classifies each as decided, provisionally decided, deferred or blocking verification. It is still a
proposal: no source files, directories or build files exist, no unit has started, and Unit 1 does not start until
the owner has decided Q6–Q9 and approved the V1–V3 results. The document must be rewritten to describe what actually exists as units land.

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
  is not analysed here. Unit 8 must decide the criterion explicitly.
- `tol_residual = 1e-10` is declared in the options (`device.py`); where the C++ solver uses it was not read.
- The legacy `Device3D` scales from whatever array it was built with (described earlier in this section); the explicit
  `Ns` input of the new design is the fix.
- The constants used by the C++ core (`core/include/tcad/physics/materials.hpp`: `kQ = 1.602176634e-19`,
  `kKB = 1.380649e-23`, `kEPS0 = 8.8541878128e-14`) are the same values as `pytcad/constants.py`.

### 6.2 Scaled Newton variables and tolerances

Newton works on scaled unknowns (ψ, then n and p in the form the legacy uses). Legacy defaults:
`tol_update = 1e-8` (maximum scaled update), `tol_residual = 1e-10`, `max_dpsi = 5.0` (damping cap on the
scaled potential update) (`pytcad/device.py`, solver options dataclass [verified]). Tolerances are stated in
scaled units. The public API stays in the cm/V convention; scaling is internal to `assemble`/`solve`.
Whether the new solver reproduces these exact defaults is to be confirmed unit by unit.

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
  during Newton overshoot must not be applied to the final converged value, and a linear solve whose relative
  residual exceeds a threshold is an error (legacy: 1e-6 for PARDISO [verified, `CLAUDE.md`]).
- `std::expected` needs `<expected>` in the pinned MSVC toolset; this is part of V1.

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
3. Backend choice and options (tolerance, iteration limit, preconditioner, thread count) are explicit
   configuration values, not environment variables. The legacy used environment variables such as
   `PYTCAD_LINSOLVE_BACKEND` (`CLAUDE.md` [verified]); they are not carried over.
4. Every solve returns `std::expected` with diagnostics (relative residual, and iteration count where it
   applies). Direct and iterative backends share one result shape. A residual above the threshold is an error even
   if the backend itself reported success (legacy: PARDISO perturbs tiny pivots instead of failing [verified]).
5. Thread count is part of backend configuration and defaults to 1 (6.8).

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

## 7. Language placement (D6: C and Fortran deferred)

| Component | Language | Reason |
|---|---|---|
| base, linalg, mesh, physics, device, results, assemble, solve, analysis, render, app | C++23 | primary language |
| Any later low-level kernel | C or Fortran | deferred; only for a measured, justified component named by the owner |

C and Fortran integration is deferred. If introduced later, it is used only for a measured and justified
component, behind a stable boundary (a C ABI), and the primary architecture stays C++23. Such a component must
be toolchain-compatible with D1 (MSVC); no Fortran compiler is chosen now. The legacy has no Fortran or C
source [verified: no `.f90/.c/.h` files in the legacy tree listing].

## 8. Proposed repository layout (nothing is created)

```text
ARCHITECTURE.md
CMakeLists.txt              top-level build; enforces C++23                 (Unit 1)
vcpkg.json                  vcpkg manifest with pinned baseline (D8)        (Unit 1)
include/NiTCAD/<layer>/...  public headers per layer                        (created per unit)
src/<layer>/...             implementations                                 (created per unit)
tests/<layer>/...           Catch2 tests mirroring src/                     (created per unit)
app/                        Win32 entry point                               (created at L8 only)
```

The layer names are decided (D5). The casing of the include root and the namespace nesting are provisional (Q1).
Not part of the structure until the owner approves them in a unit request: `CMakePresets.json`, `cmake/`
(toolchain modules), `benchmarks/` (measured performance cases) and CI configuration. D8 requires a
reproducible CI dependency and cache strategy. That strategy is a design item to specify before CI exists
(V2); it is not implemented here. A unit adds its headers, sources and tests together.

## 9. Build, test and warnings

- **Build system:** CMake. The C++23 requirement is enforced by the build, not trusted to defaults: set
  `CMAKE_CXX_STANDARD 23`, `CXX_STANDARD_REQUIRED ON`, `CXX_EXTENSIONS OFF`, and fail at configure time if the
  compiler cannot compile a probe using the C++23 features the code actually uses (for example `<expected>`).
- **Toolchain (D1, provisional):** MSVC, with a specific toolset version pinned. Verified on this machine
  (V1, 2026-10-02): Visual Studio 18 Community, MSVC toolset directory `14.51.36231`, compiler 19.51.36260
  (`_MSC_FULL_VER` 195136260). **There is no stable `/std:c++23` switch in this toolset:** `cl` lists only
  `c++14|c++17|c++20|latest`, and `/std:c++23` is accepted with warning D9002 ("ignoring unknown option") and
  silently falls back to C++14 (`_MSVC_LANG` = 201402). `/std:c++23preview` gives `_MSVC_LANG` 202302 and
  `/std:c++latest` gives 202400. CMake 4.3.1 maps `CXX_STANDARD 23` to `-std:c++latest`. Which of the two modes to
  pin is open (Q6).
- **Dependencies (D8, provisional):** CMake with vcpkg in manifest mode: a pinned vcpkg baseline, explicit
  version constraints and overrides where required, and a reproducible CI dependency and cache strategy.
  Verified locally (V2): a manifest with `builtin-baseline`, `version>=` and `overrides` resolved and built Eigen
  and Catch2 v3 with the Visual Studio-bundled vcpkg, and an Eigen SparseLU test plus a Catch2 test passed under
  `-std:c++latest`. The CI cache strategy is not specified (Q9).
- **Warning policy** (defined here, specific flags set at Unit 1 after V1): first-party code is compiled at a high
  warning level with warnings treated as errors; third-party headers are included as system headers so their
  warnings do not fail the build.
- **Test levels** per component, from the legacy practice: analytic and limiting cases → published-value
  regression → Jacobian vs finite differences → dimensional-reduction identity → convergence and mesh
  independence → benchmark. A compile-only check is never the sole gate for a numerical unit.

## 10. Validation criteria pinned from the legacy

These are the acceptance criteria the new units should meet. Source and verification status are stated for each.

| Criterion | Value | Source / status |
|---|---|---|
| Constants | Values as in `pytcad/constants.py`: q = 1.602176634e-19 C, k = 1.380649e-23 J/K, ε₀ = 8.8541878128e-14 F/cm, ħ = 1.054571817e-34 J·s, m₀ = 9.1093837015e-31 kg | [verified] file contents. The file does not name the CODATA edition; these numbers are the CODATA 2018 values per the author's knowledge [unverified: edition]. The new base unit should pin CODATA 2018 explicitly. |
| Thermal voltage | V_T(300 K) = 0.025852 V | [derived] 0.0258519998 V from the legacy k and q |
| n_i(300 K), Si | ≈ 1.0674e10 cm⁻³ | [derived] hand-computed from legacy `materials.py` (Eg0 = 1.17 eV, α = 4.73e-4 eV/K, β = 636 K, Nc300 = 2.86e19, Nv300 = 3.10e19, k = 8.617333262e-5 eV/K): Eg(300) = 1.12452 eV, n_i = 1.06738e10. **Not a legacy test assertion.** The legacy test only requires 9e9 < n_i < 1.6e10 (`tests/test_model_benchmarks.py`, `test_ni_300k_within_accepted_band` [verified]). Proposed new gate: 1.0674e10 with relative tolerance 1e-4, which holds if the same parameters are used. |
| Built-in potential | within 2e-3 V of V_T ln(Nd·Na/n_i²) | [verified] `tests/test_device1d_native_gates.py` (1e17/1e17 fixture) |
| Diode J(0.5 V) | 1.280e-2 A/cm² within 1% | [verified] `tests/test_device1d_native_gates.py::test_g2_forward_current_matches_documented_value`, fixture `_diode()` in that file. The same test records the short-base analytic value 1.321e-2. Fixture parameters are to be extracted when that unit is requested. |
| Diode law | J/J_ideal in (0.85, 1.15) at 0.5 V; ideality within ±0.02 of 1 for V ≥ 0.3 V | [verified] `tests/test_validation.py::test_ideal_diode_law`. That test uses Caughey–Thomas mobility and has SRH on by default (`Models` dataclass, `device.py`: `doping_mobility=True`, `srh=True` [verified]); this is why both come before the bias solve. |
| Jacobian vs finite differences | worst per-column-normalized error ≤ 5e-5 over ≥ 80 columns | [verified] `tests/test_m13_solver.py` ("house FD-Jacobian gate (<= 5e-5, >= 80 columns)"). The exact normalization is to be read from that file at the unit. |
| Dimensional reduction | a transversely uniform 2D solution is y-independent to atol 1e-9 | [verified] `tests/test_validation_2d.py`. Legacy also reports 1.11e-16 V for 3D → 2D (`examples/05_3d_reduces_to_2d.py`, quoted in legacy `ARCHITECTURE.md` [unverified: not re-run]). New tolerances are set per unit. |
| Mesh geometry | total volume/length matches the domain to 1e-14 (2D) / 1e-10 (3D) | [verified] `tests/test_validation_2d.py`, `tests/test_validation_3d.py` |

## 11. Proposed build order (each unit = one branch)

Unit 1 does not start until V1, V2 and V3 (section 14) are resolved. The order is by dependency, so each unit is
testable on arrival. Units 4–9 are dimension-generic from the start (D4) and use no tensor-grid assumptions (6.11).

| # | Unit | Layer | Legacy reference | Gate |
|---|---|---|---|---|
| 1 | Build scaffold: CMake, C++23 gate, vcpkg manifest with pinned baseline, Catch2 v3 harness | — | `core/CMakeLists.txt` (reference for options only) | configure fails on a toolset without the needed C++23 features; dependencies resolve from the pinned manifest; one trivial test runs |
| 2 | Constants, units, `Error` type | base | `constants.py`, `core/include/tcad/base/errors.hpp` | CODATA 2018 values; V_T(300 K) = 0.025852 V |
| 3 | Sparse matrix, backend-neutral solver interface (6.10), Eigen SparseLU backend | linalg | `linsolve.py`, `core/src/solver/direct_lu.cpp` | known systems (e.g. analytic tridiagonal); relative residual check; singular system returns an error; symbolic-reuse path exercised. A test-only second backend to prove neutrality is proposed, owner to confirm (Q5) |
| 4 | Generic node/edge/control-volume mesh; tensor-grid constructors for D = 1, 2, 3 | mesh | `mesh.py`, `mesh2d.py`, `mesh3d.py`, `core/include/tcad/mesh/stencil.hpp` | total volume matches domain (section 10); positive dual volumes; consistent edge geometry across D = 1, 2, 3 on a uniform grid; public interface contains no (i, j, k) indexing |
| 5 | Si material parameters, Caughey–Thomas mobility, SRH recombination, Boltzmann statistics (value + partials) | physics | `materials.py`, `core/include/tcad/physics/` | n_i(300 K) ≈ 1.0674e10; published mobility values (legacy `test_caughey_thomas_matches_published_silicon_values`); derivatives vs finite differences; SRH vanishes at equilibrium |
| 6 | Device description and ohmic contact data | device | `core/include/tcad/device1d/inputs.hpp` (reference only) | construction and validation; invalid input returns an error |
| 7 | Scaling (V3), SG/Bernoulli flux, equilibrium Poisson residual + Jacobian, ohmic boundary | assemble | `core/src/device1d/device1d.cpp`, `core/src/device1d/inputs.cpp` | Bernoulli limits and symmetry; FD-Jacobian gate (section 10) |
| 8 | Newton solver (scaled variables) and equilibrium solve | solve | `device.py` options, `device1d.cpp` | built-in potential within 2e-3 V; bulk neutrality; convergence; non-convergence returns an error value |
| 9 | Electron/hole continuity assembly and bias solve; first end-to-end gate | assemble, solve | `device1d.cpp`, `tests/test_validation.py`, `tests/test_device1d_native_gates.py` | J(0.5 V) = 1.280e-2 A/cm² ± 1%; ideal-diode law; current continuity; mesh independence; uniform 2D/3D reproduces 1D to a tolerance set at this unit |
| 10 | Result representation, cancellation and progress | results, solve | — (new) | cancellation returns partial results; progress is monotonic |
| 11+ | Everything else: Fermi–Dirac, Auger, band-gap narrowing, field mobility, heterojunctions, impact ionization, BTBT, transient, AC, thermal, process, unstructured meshes, PARDISO/iterative backends, file formats, analysis, render, app | deferred | per the audit | per unit, when the owner requests |

Units 1–9 are the smallest end-to-end vertical slice: a validated drift-diffusion diode. Everything past
Unit 10 is deferred until sequenced. Nothing in this table is started.

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

**Derived here, not a legacy assertion**
- V_T(300 K) = 0.0258519998 V; n_i(300 K) = 1.06738e10 cm⁻³ (inputs in section 10).

**Unverified: do not rely on these**
- That the constants file corresponds to CODATA 2018 (edition not named in the file).
- The exact normalization of the 5e-5 Jacobian gate; the unit of `D0_REF`; where `tol_residual` is used.
- clang-cl and MinGW claims about C++23 (not tested; only MSVC was probed). Behaviour of other MSVC toolsets or
  a CI runner's Visual Studio version.
- Fortran toolchain claims: that MSVC has no Fortran compiler; which Fortran compilers (Intel ifx, LLVM Flang,
  gfortran) can interoperate with which C++ ABI on Windows; CMake generator support for each.
- The MKL licence terms. Whether a CI runner can reproduce the local vcpkg build (cache behaviour).
- Whether Direct3D 12 / Direct2D headers work under a MinGW toolchain.

## 13. Risks

1. **C++23 on Windows.** Compiler mode and library support are unverified (section 12). Verified: no stable
   `/std:c++23` exists in the pinned toolset family, and an unknown `/std` value silently falls back to C++14.
   Mitigation: an explicit mode choice (Q6), a configure-time check of `_MSVC_LANG` and the feature macros, and
   treating D9002 as an error.
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

## 14. Decisions

Recorded by the owner. D7 of an earlier draft (the branch name for this file) was not a decision and is
withdrawn; numbering is otherwise kept stable. This proposal lives on `architecture/architecture-proposal`.

### 14.1 Decided

| # | Decision |
|---|---|
| D2 | Catch2 v3 is the C++ test framework. Its availability through vcpkg is checked under V2; the decision stands. |
| D4 | Dimension-generic architecture for 1D, 2D and 3D. |
| D5 | Root namespace `NiTCAD`. Layers: `base`, `linalg`, `mesh`, `physics`, `device`, `assemble`, `solve`, `results`, `analysis`, `render`, `app`. (Namespace nesting and include-root casing: provisional, Q1.) |
| R1 | `std::expected` for recoverable errors; exceptions only for exceptional or programmer-error cases; no exceptions in hot numerical kernels or Newton iterations (6.7). |
| R4 | Boltzmann carrier statistics initially; Fermi–Dirac and other statistics deferred until requested (section 5, item 7). |

### 14.2 Provisionally decided (stand unless verification or Unit work disproves them)

| # | Decision | Provisional because |
|---|---|---|
| D1 | MSVC is the Windows C++23 toolchain, with one specific toolset pinned. C++ compiler only. | V1 is done, but it found no stable C++23 mode in this toolset (Q6, Q9). |
| D3 | Eigen SparseLU is the initial sparse direct solver behind a backend-neutral interface that allows PARDISO and iterative backends later without redesign (6.10). | The exact interface shape is settled at Unit 3. |
| D8 | CMake + vcpkg manifest mode, pinned baseline, explicit version constraints/overrides where required, reproducible CI dependency/cache strategy. | V2 is done for a local build; the baseline choice (Q7), triplet (Q8) and CI cache strategy (Q9) are open. |
| D5 (part) | `NiTCAD::<layer>` nesting; `include/NiTCAD/<layer>/`. | Owner named the root namespace and layers only (Q1). |

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
| V1 | Toolset: VS 18 Community, MSVC dir `14.51.36231`, compiler 19.51.36260. Required features are present (`std::expected`, `std::mdspan`, explicit object parameters, multidimensional subscript, `std::println`, `<stop_token>`/`jthread`) under `/std:c++23preview` and `/std:c++latest`. **No stable `/std:c++23`**; it is silently ignored (D9002, falls back to C++14). CMake 4.3.1 emits `-std:c++latest` for `CXX_STANDARD 23`. Details in section 12. | Verified. **Finding needs an owner decision (Q6).** |
| V2 | The Visual Studio-bundled vcpkg works in manifest mode with a pinned baseline, `version>=` constraints and overrides. Eigen (5.0.1 at the test baseline, 3.4.0 by override) and Catch2 3.16.0 built and passed a probe test under `-std:c++latest`. Not verified: CI reproduction and cache behaviour, a standalone vcpkg checkout, whether the redirect of vcpkg's build/download/package directories out of Program Files was necessary (done as a precaution). | Verified locally. **Open items: Q7, Q8, Q9.** |
| V3 | Scaling definitions, scaled variables and Newton tolerances read from `inputs.cpp` and `device1d.cpp` and recorded in 6.1. One open question remains for Unit 8 (convergence criterion on the clipped correction). | Verified |

Verifications due at their own unit, not blocking Unit 1: CODATA edition of the constants (Unit 2), the exact
normalization of the 5e-5 Jacobian gate (Unit 7), the `_diode()` fixture parameters (Unit 9).

### 14.5 Remaining architectural questions

| # | Question | Needed by |
|---|---|---|
| Q1 | Namespace nesting (`NiTCAD::<layer>` vs one flat namespace) and include-root casing (`include/NiTCAD/` vs lowercase) | Unit 1 |
| Q2 | Mechanism for programmer-error preconditions: exception, assertion or terminate (6.7) | Unit 2 |
| Q3 | Whether to use the Visual Studio-bundled vcpkg (verified working, V2) or a pinned standalone checkout (untested), and whether `CMakePresets.json` and a CI workflow are part of Unit 1 | Unit 1 |
| Q4 | Final `Error` category list (6.7) | Unit 2 |
| Q5 | Whether Unit 3 includes a test-only second backend to prove backend neutrality | Unit 3 |
| Q6 | Which MSVC mode enforces C++23: `/std:c++23preview` (`_MSVC_LANG` 202302) or `/std:c++latest` (202400, rolling, what CMake emits by default). Either way the build should set the flag explicitly and fail configuration if `_MSVC_LANG` or the feature macros do not match. Candidate: `/std:c++23preview`, revisited when a stable `/std:c++23` exists | Unit 1 |
| Q7 | vcpkg baseline commit and update policy; Eigen major version (the baseline tested gives 5.0.1; the legacy Eigen version was not checked, and an override to 3.4.0 also worked) | Unit 1 |
| Q8 | vcpkg triplet (the probe used `x64-windows-static-md`, chosen by me for the test only) | Unit 1 |
| Q9 | How the toolset is pinned (VS version, toolset `14.51.36231`, edition and licence of Visual Studio Community, how CI provisions the same toolset) and the CI vcpkg binary-cache strategy | Unit 1 |

**Approval workflow:** the owner names each unit. No unit starts before Q6–Q9 are decided and the owner
approves the V1–V3 results.
