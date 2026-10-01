# NT-SemTCAD — Architecture (PROPOSAL, awaiting decisions)

**Status:** proposal only. No source files, directories or build files exist yet. Nothing below is
final until the owner approves it in section 9. This document describes intent; it must be updated
to describe what actually exists as each unit lands.

Reference: the legacy repository `NT-SemTCAD` (science, algorithms, tests and benchmarks only; its
architecture is not carried over).

---

## 1. Fixed constraints (already decided)

| Item | Decision |
|---|---|
| Platform | Windows only |
| Languages | C++23 (primary), C and Fortran only where technically justified |
| C++ standard | C++23, enforced by the build system rather than compiler defaults |
| Desktop layer | Win32 + Direct3D 12 + Direct2D. No Qt, VTK, Python GUI or other GUI framework |
| Runtime | No Python anywhere in the runtime |
| Tests | Written in C/C++ |
| Git | One branch per architectural unit, no auto-merge, `main` stays clean |
| Scope | Feature freeze: only what is explicitly requested is built |

## 2. Design principles

1. One responsibility per file; explicit dependencies; no global state; no circular dependencies.
2. Dependencies point one way only (section 3). A lower layer never includes a higher one.
3. **Physics is separate from numerics, and both are separate from the device.** The legacy
   `device.py` / `device2d.py` / `device3d.py` each mixed geometry, scaling, physics, assembly,
   Jacobian and solve. The new design splits those apart.
4. **Dimension is a parameter, not a code path.** Legacy 1D/2D/3D were separate code with
   dimensional-reduction tests holding them together. The proposal is one dimension-generic
   discretization, with the legacy reduction identities (3D uniform → 2D → 1D) kept as test gates.
5. **The solver core has no UI knowledge and the UI has no physics knowledge.** They meet only
   through plain data (mesh, fields, results).
6. Scientific behaviour is preserved: equations, units (cm, cm⁻³, V, A/cm²), constants, sign
   conventions, tolerances. Architecture changes do not change physics.

## 3. Proposed layers (dependency order, lowest first)

```text
 L7  app          Win32 application: windows, input, document, run control
 L6  render       Direct3D 12 / Direct2D drawing of meshes, fields, plots
 L5  analysis     observables derived from results (I-V, C-V, Vth, ...)
 L4  solve        nonlinear/continuation drivers, sweeps, run results
 L3  assemble     residual + Jacobian assembly from physics over a mesh
 L2  physics      models: statistics, mobility, recombination, ionization, ...
 L1  mesh/numerics  mesh data, discretization (Scharfetter-Gummel), linear algebra
 L0  base         constants, units, errors, checked numeric utilities
```

Allowed dependencies: each layer may use only layers below it. `render` and `app` may use `analysis`
and `solve` results but never include physics or assembly headers. `analysis` consumes results only.

Rationale for each cut:

| Layer | Owns | Does NOT own |
|---|---|---|
| base | physical constants, unit convention, error type | any physics |
| mesh/numerics | geometry and mesh, SG flux, sparse linear solve | device meaning |
| physics | pure functions of local state (n, p, T, field, doping) with cited parameters | meshes, matrices |
| assemble | turning physics plus mesh into F(x) and J(x) | Newton iteration |
| solve | Newton, damping, continuation, bias sweeps | UI, plotting |
| analysis | derived quantities from stored results | solving |
| render / app | pixels, windows, user interaction | numerics |

## 4. Proposed language placement

| Component | Language | Reason |
|---|---|---|
| base, mesh, assemble, solve, analysis, render, app | C++23 | primary language; owns architecture |
| Stable low-level numeric kernels (if profiling later justifies) | C or Fortran | only after measurement; not assumed |
| Dense/sparse linear solve | **decision D3** | a vendored solver vs a Windows system library |

Recommendation: start everything in C++23. Introduce Fortran or C only for a specific kernel, when
you name it and a benchmark justifies it. Legacy measurements show that the direct-solver
algorithm, not the implementation language, dominated 3D wall time (98% of a 24³ solve was in the
sparse LU), so the linear solver choice matters more than the assembly language.

## 5. Proposed repository layout (NOT created; for your approval)

```text
ARCHITECTURE.md
CMakeLists.txt                  top-level build; enforces C++23 and the warning policy
cmake/                          toolchain notes/modules, only if needed
include/tcad/<layer>/...        public headers per layer
src/<layer>/...                 implementations
tests/<layer>/...               C++ tests, mirroring src/
benchmarks/                     measured performance cases (added when first needed)
app/                            Win32 entry point (added at L7)
```

The layer directory names (`base`, `mesh`, `physics`, `assemble`, `solve`, `analysis`, `render`,
`app`) are suggestions. You said you control the tree, so no directory is created until you approve
its name. A unit is created as headers, source and tests together, on its own branch.

## 6. Build and test (proposal)

- **Build system:** CMake with `CMAKE_CXX_STANDARD 23`, `CXX_STANDARD_REQUIRED ON`, extensions off.
  A configure-time check fails if the compiler cannot do C++23. This is the mechanism that enforces
  the C++23 requirement.
- **Compiler:** **decision D1.**
- **Test framework:** **decision D2.**
- **Test levels**, selected per component (from the legacy repo's practice):
  analytic/limiting cases → published-value regression → Jacobian vs finite differences (tolerance
  5e-5 in the legacy suite) → dimensional-reduction identity → convergence/mesh independence →
  benchmark. Never "it compiled" as the only gate for numerics.
- Bit-identity digests are **not** carried over. The legacy repo found they pin one machine's
  floating-point summation order. Use tolerances instead.

## 7. Proposed build order (each item = one branch)

The order follows the dependency graph so every unit is testable on arrival. Items marked
*[optional]* are legacy capabilities that stay frozen out until you request them.

| # | Unit | Layer | Legacy reference | First validation |
|---|---|---|---|---|
| 1 | Build scaffold (CMake, C++23 gate, test harness) | — | none (new) | configure fails on a non-C++23 compiler |
| 2 | Constants, units, errors | base | `constants.py` | CODATA values, V_T = 25.852 mV at 300 K |
| 3 | Material parameters (Si first) | physics | `materials.py` | n_i(300 K) in accepted band |
| 4 | 1D mesh, SG flux (Bernoulli) | mesh | `mesh.py`, `device.py` | Bernoulli limits and symmetry |
| 5 | Equilibrium Poisson (1D) | assemble/solve | `device.py` | built-in potential, neutrality, FD-Jacobian |
| 6 | Drift-diffusion bias solve (1D) | assemble/solve | `device.py` | ideal-diode law, current continuity |
| 7 | Mobility, SRH/Auger, BGN | physics | `materials.py` | published Caughey-Thomas/Slotboom values |
| 8 | Dimension-generic mesh and 2D | mesh | `mesh2d.py`, `device2d.py` | 2D uniform reproduces 1D |
| 9 | 3D | mesh | `device3d.py` | 3D uniform reproduces 2D |
| 10 | Linear solver strategy | mesh/numerics | `linsolve.py` | residual check; scaling benchmark |
| 11+ | Everything else (heterojunctions, Fermi–Dirac, impact ionization, BTBT, transient, AC, thermal, process, …) | *[optional]* | per audit | per unit |
| L6/L7 | Win32 / D3D12 / D2D application | render, app | none | per unit |

Items 1–6 form the smallest end-to-end vertical slice: a validated 1D diode. Everything past item 6
is deferred until you sequence it.

## 8. Risks

1. **Linear solver on Windows.** The legacy repo used MKL PARDISO (loaded at runtime) and Eigen
   SparseLU, with optional PETSc. The new solver choice affects licensing and the build.
2. **Fortran toolchain.** The legacy repo contains no Fortran. If it is wanted, a Fortran compiler
   and C++/Fortran interop must be set up and tested on Windows.
3. **Dimension-generic design cost.** One generic discretization is cleaner but is a larger initial
   design effort than three direct ports. The legacy 1D/2D/3D agreement is the safety net.
4. **Unvalidated legacy physics.** Several legacy models are disclosed simplifications or lack a
   published reference (hydrodynamic closure, MC implant, TED, self-heating). They are deferred, and
   will be reimplemented only with an explicit scope statement.
5. **Direct3D 12 verbosity.** The render layer is large and low-level and does not affect solver
   correctness. It is deliberately last.

## 9. Decisions needed from the owner

| # | Question | Options | Recommendation |
|---|---|---|---|
| D1 | Compiler/toolchain | MSVC (VS Build Tools) · clang-cl · MinGW-w64 g++ | MSVC: native to Win32/D3D12/D2D and best Windows tooling. Legacy used MinGW, so this is a change, not a continuation. |
| D2 | C++ test framework | Catch2 · doctest · GoogleTest · plain `ctest` + asserts | Catch2 (the legacy core already planned it) |
| D3 | Linear solver | Eigen (header-only) · MKL PARDISO · hand-written | Eigen first. Add PARDISO only when a benchmark justifies it. |
| D4 | Dimension-generic mesh vs per-dimension | generic · separate 1D/2D/3D | generic, with reduction-identity tests |
| D5 | Layer and directory names (section 5) | as proposed · yours | your naming |
| D6 | Fortran/C | defer entirely · name a first kernel now | defer until a measured need |
| D7 | Branch for this file | `architecture/architecture-proposal` (used) | — |

Please answer D1–D6 (or edit this file directly), and then name the first unit. Unit 1 (build
scaffold) is the natural first request, but I will not create it until you say so.
