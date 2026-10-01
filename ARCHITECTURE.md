# NiTCAD — Architecture

**PROPOSAL, awaiting decisions.**

This document is a proposal. No source files, directories or build files exist yet. Nothing in it is
final until the owner approves the architecture and decisions D1–D6 and D8 (section 14). It must be
rewritten to describe what actually exists as units land.

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

## 4. Layer responsibilities

| Layer | Owns | Does NOT own |
|---|---|---|
| base | physical constants, unit convention, error types | physics |
| linalg | sparse matrix storage, baseline direct solve (D3), residual check | meshes, Newton |
| mesh | node/edge/control-volume graph, per-edge length and coupling area, per-node volume; tensor-grid constructors for 1D/2D/3D | doping, materials, linear algebra |
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

## 6. Cross-cutting design

### 6.1 de Mari scaling

The legacy solves in scaled variables. Its `Inputs` struct carries `VT, ni, Ns, LD, J0, R0, eps`, with a
doping scale `Ns` that defaults to `max(|doping|, ni)` unless overridden
(`inputs.hpp`; `Device1D` comment on `Ns` [verified]). The legacy `Device3D` derives its scaling from the
maximum |doping| of whatever array it was built with, so two devices covering different slices silently
disagreed on units until an `Ns_override` was added (`CLAUDE.md`, "Physics/model conventions" [verified]).

Proposal: scaling is an explicit `assemble` input derived once per problem, never recomputed per slice
or hidden inside a device. Exact definitions of `LD`, `J0` and `R0` are not asserted here. They will be read
from `core/src/device1d/device1d.cpp` when the scaling unit is requested.

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
on-disk input format is **not decided** (remaining decision R2). The legacy wire format was a Python JSON
`DeviceSpec` and is not carried over.

### 6.6 Result representation

`results` holds fields indexed by mesh node, terminal quantities per bias point, and a run record
(inputs hash/identity, options, per-iteration convergence history, converged flag). It is plain data.
An on-disk result format is **not decided** (R2).

### 6.7 Error policy

Proposal: one `tcad::Error` base type with typed subclasses (degenerate mesh, singular linear system,
non-convergence, invalid input), failing loudly on invalid input instead of silently clamping. The
mechanism (exceptions vs `std::expected`) is **not decided** (R1). The legacy uses a typed C++
exception hierarchy in `core/include/tcad/base/errors.hpp` [verified]. Two legacy rules are kept: a
clamp used during Newton overshoot must not be applied to the final converged value, and a linear solve whose relative
residual exceeds a threshold is an error (legacy: 1e-6 for PARDISO [verified, `CLAUDE.md`]).

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

## 7. Language placement

| Component | Language | Reason |
|---|---|---|
| base, linalg, mesh, physics, device, results, assemble, solve, analysis, render, app | C++23 | primary language |
| Any later low-level kernel | C or Fortran | only when the owner names the kernel and a benchmark justifies it |

Starting everything in C++23 is the recommendation (D6). The legacy has no Fortran or C source [verified:
no `.f90/.c/.h` files in the legacy tree listing]. A language is not added merely to use all three.

## 8. Proposed repository layout (nothing is created)

```text
ARCHITECTURE.md
CMakeLists.txt        top-level build; enforces C++23              (Unit 1)
include/<ns>/<layer>/...   public headers per layer                (created per unit)
src/<layer>/...            implementations                         (created per unit)
tests/<layer>/...          C++ tests mirroring src/                (created per unit)
app/                       Win32 entry point                       (created at L8 only)
```

Not part of the proposed structure until actually required: `cmake/` (toolchain modules), `benchmarks/`
(measured performance cases), CI configuration. The layer directory names and the namespace (`tcad` in the
legacy) are suggestions pending D5. A unit adds its headers, sources and tests together.

## 9. Build, test and warnings

- **Build system:** CMake. The C++23 requirement is enforced by the build, not trusted to defaults:
  set `CMAKE_CXX_STANDARD 23`, `CXX_STANDARD_REQUIRED ON`, `CXX_EXTENSIONS OFF`, and fail at configure time
  if the compiler cannot compile a probe using the C++23 features the code actually uses. Whether the
  chosen compiler's "C++23" mode is a stable mode or a preview/latest mode, and which library features it
  supports, **requires verification** before D1 is approved (section 12).
- **Warning policy** (defined here, specific flags depend on D1): first-party code is compiled at a high
  warning level with warnings treated as errors; third-party headers are included as system headers so
  their warnings do not fail the build.
- **Compiler, test framework, solver, dependency mechanism:** see D1, D2, D3, D8.
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

D3 and D8 (and D1/D2, which the scaffold needs) must be approved before Unit 1 starts or, at the latest, before
the unit that uses them. The order is by dependency, so each unit is testable on arrival. Items 4–9 are
dimension-generic from the start (D4).

| # | Unit | Layer | Legacy reference | Gate |
|---|---|---|---|---|
| 1 | Build scaffold: CMake, C++23 gate, test harness, dependency mechanism (D8) | — | `core/CMakeLists.txt` (reference for options only) | configure fails on a compiler without the needed C++23 features; one trivial test runs |
| 2 | Constants, units, error base | base | `constants.py`, `core/include/tcad/base/errors.hpp` | CODATA 2018 values; V_T(300 K) = 0.025852 V |
| 3 | Sparse matrix + baseline direct solve (D3) | linalg | `linsolve.py`, `core/src/solver/direct_lu.cpp` | known systems (e.g. analytic tridiagonal), relative residual check, singular system raises |
| 4 | Generic node/edge/control-volume mesh; tensor-grid constructors for D = 1, 2, 3 | mesh | `mesh.py`, `mesh2d.py`, `mesh3d.py`, `core/include/tcad/mesh/stencil.hpp` | total volume matches domain (section 10); positive dual volumes; 1D/2D/3D give consistent edge geometry on a uniform grid |
| 5 | Si material parameters, Caughey–Thomas mobility, SRH recombination (value + partials) | physics | `materials.py`, `core/include/tcad/physics/` | n_i(300 K) ≈ 1.0674e10; published mobility values (legacy `test_caughey_thomas_matches_published_silicon_values`); derivatives vs finite differences; SRH vanishes at equilibrium |
| 6 | Device description and ohmic contact data | device | `core/include/tcad/device1d/inputs.hpp` (reference only) | construction and validation; invalid input raises |
| 7 | Scaling, SG/Bernoulli flux, equilibrium Poisson residual + Jacobian, ohmic boundary | assemble | `core/src/device1d/device1d.cpp`, `core/src/device1d/inputs.cpp` | Bernoulli limits and symmetry; FD-Jacobian gate (section 10) |
| 8 | Newton solver (scaled variables) and equilibrium solve | solve | `device.py` options, `device1d.cpp` | built-in potential within 2e-3 V; bulk neutrality; convergence |
| 9 | Electron/hole continuity assembly and bias solve; first end-to-end gate | assemble, solve | `device1d.cpp`, `tests/test_validation.py`, `tests/test_device1d_native_gates.py` | J(0.5 V) = 1.280e-2 A/cm² ± 1%; ideal-diode law; current continuity; mesh independence; uniform 2D/3D reproduces 1D to a tolerance set at this unit |
| 10 | Result representation, cancellation and progress | results, solve | — (new) | cancellation returns partial results; progress is monotonic |
| 11+ | Everything else: Fermi–Dirac, Auger, band-gap narrowing, field mobility, heterojunctions, impact ionization, BTBT, transient, AC, thermal, process, analysis, render, app | deferred | per the audit | per unit, when the owner requests |

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

**Derived here, not a legacy assertion**
- V_T(300 K) = 0.0258519998 V; n_i(300 K) = 1.06738e10 cm⁻³ (inputs in section 10).

**Unverified: do not rely on these**
- That the constants file corresponds to CODATA 2018 (edition not named in the file).
- The exact definitions of `LD`, `J0`, `R0`, and the exact normalization of the 5e-5 Jacobian gate.
- Any MSVC/clang/MinGW claim about C++23 mode, library feature support (`std::expected`, `std::mdspan`,
  explicit object parameters), multidimensional subscript support, or Eigen compatibility with the C++23 flag.
- Fortran toolchain claims: that MSVC has no Fortran compiler; which Fortran compilers (Intel ifx, LLVM Flang,
  gfortran) can interoperate with which C++ ABI on Windows; CMake generator support for each.
- Catch2 v3 packaging (compiled library vs header-only), and the MKL licence terms.
- Whether Direct3D 12 / Direct2D headers work under a MinGW toolchain.

## 13. Risks

1. **C++23 on Windows.** Compiler mode and library support are unverified (section 12). Mitigation: verify
   before D1 is approved; the configure-time probe fails loudly.
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

D7 of the earlier draft (the branch name for this file) was not a decision and is withdrawn. It is a status
item: this proposal lives on `architecture/architecture-proposal`. Numbering is otherwise kept stable.

| # | Decision | Options | Recommendation | Status |
|---|---|---|---|---|
| D1 | **Windows C++23 toolchain**: C++ compiler and version, standard library, build generator, and the C++23 gate (the configure-time probe). The C++ compiler only; Fortran is D6 | A. MSVC (VS Build Tools, minimum version to be verified) + Ninja or VS generator · B. clang-cl on the MSVC STL · C. MinGW-w64 g++ (the legacy toolchain) | A, **subject to verifying** the C++23 mode and the library features the code uses. B is a fallback if MSVC lacks one. C is not recommended: native Win32/D3D12/D2D headers under MinGW are unverified. | Open |
| D2 | C++ test framework | Catch2 · doctest · GoogleTest · plain `ctest` with asserts | Catch2 (the legacy core planned it behind `TCAD_BUILD_TESTS`); verify its packaging under D8 | Open |
| D3 | Baseline sparse direct linear solver (needed before the first Newton unit, which is Unit 8; Unit 3 delivers it) | Eigen SparseLU · MKL PARDISO · own implementation | Eigen SparseLU behind a small `linalg` interface, so PARDISO can be added later if a new-repo benchmark justifies it. Not decided for the owner. | Open, **must be approved before Unit 3** |
| D4 | Dimensional strategy | A. dimension-generic node/edge graph (box method) from the first unit · B. per-dimension codebases | A. This document is written assuming A. Choosing B would change units 4–9. | Open |
| D5 | Layer directory names and the C++ namespace | as in section 3 · owner's names | owner's naming | Open |
| D6 | Fortran and C | defer entirely · name a first kernel now | Defer until a measured need exists. If adopted later: a Fortran compiler that is ABI-compatible with the D1 compiler, behind a C ABI. | Open |
| D7 | withdrawn (see above) | — | — | Withdrawn |
| D8 | **Dependency acquisition and management** (Catch2, Eigen, any later dependency) | A. vendored in the repo · B. a package manager (e.g. vcpkg, Conan) · C. CMake `FetchContent` · D. preinstalled/system dependency · E. another controlled mechanism | Not chosen here. Each option affects the owner-controlled file tree (A adds a third-party directory), reproducibility (D depends on the machine) and network use at configure time (B, C). The mechanism must pin versions. | Open, **must be approved before Unit 1** |

### Remaining architectural decisions (not part of D1–D8)

| # | Question |
|---|---|
| R1 | Error mechanism: exceptions, `std::expected`, or a mix (6.7) |
| R2 | On-disk input and result formats (6.5, 6.6); none exist yet |
| R3 | Whether unstructured meshes are in the initial scope or the mesh stays tensor-grid until requested |
| R4 | Statistics scope of the first physics unit: Boltzmann only (as in the proposed Unit 5) or also Fermi–Dirac |
| R5 | Exact scaling definitions to adopt (6.1), confirmed against `device1d.cpp` when requested |

**Approval workflow:** this document stays marked "PROPOSAL, awaiting decisions" until the architecture and
D1–D6 and D8 are approved. No unit starts before that, and the owner names each unit.
