# NiTCAD — Decisions awaiting review

**Status: Q1–Q5, Q10 and Q11 ACCEPTED by the owner on 2026-10-02** ("I like opus decision go ahead and start") and
recorded in `ARCHITECTURE.md` section 14.1. The new questions N1–N7 and the inconsistency list below were **not**
explicitly ruled on: N1–N4 are applied as provisional Unit 1 defaults (14.2), N5–N7 stay open (14.5).

Source: an Opus subagent (medium effort, read-only) was given `ARCHITECTURE.md` and the questions on 2026-10-02.
I (Claude Sonnet 5.5) have not independently verified its factual claims beyond what is marked. Claims that came from
secondary web sources, or that nobody ran, are marked **[unverified]**.

Already decided and not repeated here: D1–D6, D8, R1, R4, Q6–Q9 (see `ARCHITECTURE.md` 14.1–14.3).

## Summary

| # | Question | Recommendation | Needed before |
|---|---|---|---|
| Q1 | Namespace nesting and include casing | `NiTCAD::<layer>`; `include/NiTCAD/<layer>/` | Unit 1 |
| Q2 | Programmer-error preconditions | Always-on `NITCAD_EXPECTS` macro, fail-fast, never throws | Unit 2 |
| Q3 | vcpkg source; presets and CI in Unit 1 | Pinned standalone vcpkg; `CMakePresets.json` in Unit 1; CI later | Unit 1 |
| Q4 | `Error` categories and shape | `enum class ErrorCode` (7 values) + `struct Error` with optional context | Unit 2 |
| Q5 | Test-only second backend in Unit 3 | No; prove neutrality with a header-boundary build check and interface-only tests | Unit 3 |
| Q10 | CI toolset and cache design | GitHub-hosted Windows runner, `vcvars_ver` selection, assert exact versions, `files` binary cache | CI unit |
| Q11 | Move to an LTS toolset | Trigger on 14.52 LTS release, migrate within 60 days, re-run V1/V4 | 14.52 release |

---

## Q1 — Namespace nesting and include-root casing

- **Recommendation:** Nest each layer as `NiTCAD::<layer>`. Use `include/NiTCAD/<layer>/`, spelled exactly like the namespace.
- **Why:**
  - Nesting prevents name clashes and ADL surprises between layers (for example `mesh::Node` vs a future `render::Node`; `linalg::Matrix` vs an Eigen or D3D matrix).
  - Nesting lets a later check enforce the layer table (section 3) by namespace.
  - Windows ignores case, so the risk is consistency, not function. One spelling avoids an include that works locally and breaks on a case-sensitive checkout or cache.
  - Mixed case matches the decided root namespace.
- **Rejected:** a flat `NiTCAD` namespace (loses the layer boundary); lowercase `include/nitcad/` (a second spelling with no benefit on Windows).
- **Risks:** Git `core.ignorecase=true` can hide case-only renames, so settle the casing before any file exists.
- **Needed before:** Unit 1.
- **Owner decision:** [x] accept (2026-10-02)  [ ] change  [ ] defer

## Q2 — Enforcing programmer-error preconditions

- **Recommendation:** One project macro, `NITCAD_EXPECTS(cond)`. Always on in debug and release. On violation it logs file, line and condition, then calls `std::abort` or `__fastfail`. It never throws.
- **Why:**
  - Fits R1: no exceptions in hot kernels. The check costs one predictable branch.
  - A violated precondition is a bug, not a state to recover from; throwing invites catch-and-continue.
  - Always on keeps loud failure (6.7). A separate debug-only assert can be added later for expensive invariants.
  - Catch2 has no death tests, so preconditions are tested by testing the predicate functions directly. A test-only option routing the handler to a throwing hook is possible; keep it out of shipped builds.
- **Rejected:** exceptions (contradict R1, blur bug vs error); `assert` (disabled by `NDEBUG`); bare `std::terminate` (loses the diagnostic).
- **Risks:** C++26 contracts may later replace the macro, so keep it thin. `__fastfail` vs `abort` under a debugger **[unverified]**: not tested.
- **Needed before:** Unit 2.
- **Owner decision:** [x] accept (2026-10-02)  [ ] change  [ ] defer

## Q3 — vcpkg source; `CMakePresets.json` and CI in Unit 1

- **Recommendation:** Use a standalone vcpkg tool checkout pinned to a recorded tag or commit. Put `CMakePresets.json` in Unit 1. Keep the CI workflow out of Unit 1; add it as its own owner-named unit after Q10 is designed.
- **Why:**
  - The manifest baseline pins port versions, not the vcpkg tool. The Visual Studio-bundled tool changes with every Visual Studio update, and CI cannot assume the same copy.
  - A preset is the one place that sets generator, toolchain file, triplet, `/std:c++23preview` and `/options:strict`. Without it every developer and CI job retypes them. This is not speculative.
  - A CI workflow before the cache and toolset design exists would be speculative.
- **Rejected:** bundled vcpkg only (drifts, not reproducible on CI); vcpkg as a git submodule (adds a full vcpkg tree to the owner-controlled file tree).
- **Risks:** the standalone checkout is **[unverified]**. V2 only tested the Visual Studio-bundled one, so V2 must be re-run with a standalone checkout before Unit 1.
- **Needed before:** Unit 1 (vcpkg choice, presets). CI later.
- **Owner decision:** [x] accept (2026-10-02)  [ ] change  [ ] defer

## Q4 — `Error` categories and shape

- **Recommendation:**
  ```text
  enum class ErrorCode : std::uint8_t {
      invalid_input, degenerate_mesh, singular_system,
      inaccurate_solve, non_convergence, cancelled, resource_exhausted };
  struct Error { ErrorCode code; std::string message; std::optional<ErrorContext> context; };
  ```
  `ErrorContext` is small and plain: an optional index (node, edge, row or bias step) and an optional numeric value (residual or tolerance).
- **Why:**
  - Callers branch on the code; people read the message.
  - Typed numeric context keeps the diagnostics required by 6.10 item 4 without formatting strings inside kernels.
  - `resource_exhausted` covers an out-of-memory factorization that a backend converts at its boundary (6.7).
  - No "internal" category: internal bugs go through Q2.
- **Rejected:** `std::error_code` and categories (more machinery, string-centric); splitting `invalid_input` per layer (premature).
- **Risks:** `std::string` allocates on the error path. Acceptable because it never runs on the success path, but say so in the document.
- **Needed before:** Unit 2.
- **Owner decision:** [x] accept (2026-10-02)  [ ] change  [ ] defer

## Q5 — Test-only second backend in Unit 3

- **Recommendation:** Do not add a second backend. Prove backend neutrality with (a) a CMake test that compiles the public `linalg` headers in a translation unit with **no Eigen include path**, so any leak of an Eigen type into the interface fails the build, and (b) behavioural tests written only against the interface.
- **Why:**
  - It checks the boundary mechanically.
  - Interface-only tests run unchanged against any future backend.
  - A test-only dense or mock solver would be a speculative file the feature freeze forbids, and it proves little about PARDISO-style symbolic reuse.
- **Rejected:** a test-only dense or mock backend.
- **Risks:** whether one interface fits both direct and iterative backends stays unproven until one is requested.
- **Needed before:** Unit 3.
- **Owner decision:** [x] accept (2026-10-02)  [ ] change  [ ] defer

## Q10 — CI toolset provisioning and vcpkg binary cache

- **Recommendation:** GitHub-hosted Windows runner. Select the toolset with `vcvarsall.bat x64 -vcvars_ver=14.51`. Fail the job unless `cl` reports 19.51.36260, the SDK is 10.0.26100.0, CMake is 4.3.1 and Ninja is 1.13.2. Cache vcpkg binaries with the `files` binary source pointed at a directory saved by `actions/cache`, keyed on baseline + triplet + toolset version + `vcpkg.json` hash.
- **Why:**
  - Exact image contents cannot be pinned on hosted runners; asserting versions turns drift into a loud failure.
  - If the exact toolset is missing, the fallbacks are a self-hosted runner or a step that installs the component with the Visual Studio Installer CLI (slow).
  - The `files` source works with any cache backend.
- **Rejected:** vcpkg `x-gha` (believed removed after GitHub's cache-service change, **[unverified]**); NuGet on GitHub Packages (needs tokens and a feed; revisit if cache size limits bite).
- **Risks / unverified:**
  - Secondary sources say GitHub moved `windows-latest`/`windows-2025` to Visual Studio 2026 in June 2026 ([scivision](https://www.scivision.dev/github-actions-vs-2026), [byteiota](https://byteiota.com/github-actions-windows-vs2026-migration/)). **[unverified]**
  - **Not checked:** whether that image contains MSVC 14.51.36231 exactly, SDK 26100, CMake 4.3.1 or Ninja 1.13.2. The bundled CMake probably differs, so CMake and Ninja would be pinned separately (checksummed download or a pinned action).
  - `-vcvars_ver=14.51` selection on a runner has not been run. Hosted images update weekly, so the patch level may move.
- **Needed before:** the CI unit (not Unit 1).
- **Owner decision:** [x] accept (2026-10-02)  [ ] change  [ ] defer

## Q11 — Moving to an LTS toolset

- **Recommendation:** Trigger on general availability of MSVC 14.52 LTS (expected November 2026). Migrate within 60 days, and well before the derived end of 14.51 support (about February 2027).
- **Why:**
  - The project needs years of stability. 14.50 LTS was never tested here.
  - If 14.52 ships a stable `/std:c++23`, switch from `/std:c++23preview` at the same time.
- **Re-verification steps (one dedicated branch):**
  1. Re-run V1 (feature macros and probe program).
  2. Re-run V4 (configure gate, D9002 / `/options:strict` behaviour, Eigen 5.0.1 and Catch2 3.16.0 from the pinned baseline).
  3. Run the full test suite.
  4. Update the CI version assertions and cache key.
  5. Record the result in `ARCHITECTURE.md` section 12.
- **Rejected:** moving to 14.50 now (untested); staying on 14.51 past end of support.
- **Risks:** the 14.52 release date and its C++23 completeness are **[unverified]**. `_MSVC_LANG` under a stable flag may differ from 202302L, so the gate must accept the new value deliberately.
- **Needed before:** 14.52 release, and before 14.51 end of support.
- **Owner decision:** [x] accept (2026-10-02)  [ ] change  [ ] defer

---

## Unit 23 — Electrothermal coupling (decided by the owner on 2026-10-09)

Branch `physics/electrothermal` from `main` (`c8fd145`). The owner approved the revised scope after a self-review of the
first proposal. Legacy reference: `thermal.py`, `thermal_grid.py`, `thermal2d.py`, `thermal3d.py`,
`core/src/thermal/grid.cpp`, `tests/test_m19_thermal.py`.

| # | Decision | Owner ruling |
|---|---|---|
| T1 | Coupling | Lattice temperature is a fourth per-node Newton unknown with an exact Jacobian (monolithic). OLD: outer Gummel loop over a uniform peak-T device. NEW: per-node T. REASON: a uniform T misses local heating, and the outer loop hides thermal runaway where it fails to converge. No Gummel loop is offered. |
| T2 | Carrier flux | Full thermodynamic model (Wachutka): `J_n = -q mu_n n (grad phi_n + P_n grad T)`, likewise for holes. Thermopower with a per-material, per-carrier scattering exponent r (default -1/2, cited). With Boltzmann statistics and Nc ∝ T^1.5 the flux is linear in n (`J_n = mu_n [k T grad n + k (1+r) n grad T + n grad E_c]`), discretized by a Scharfetter–Gummel flux generalized to edge-varying T. |
| T3 | Fermi–Dirac | Included. F_0 and F_1 (r = -1/2) implemented with 40-digit validation; the Unit 14 ν-factor SG extended to edge-varying T and the thermopower. |
| T4 | Heat source | Edge divergence of the energy flux `(phi + T P) J`, minus the carrier-storage terms `(phi_n + T P_n) q dn/dt` and the hole term, so the transient and AC heat is Wachutka's local form; Joule, recombination, Peltier and Thomson heat all included. Steady energy balance exact: heat out of the thermal contacts equals sum I V to round-off. Transient: sum I V = heat + d/dt(electrostatic + carrier energy). |
| T5 | Contact Peltier heat | Metal thermopower 0; the Peltier heat `-T_c P_c I_c` is placed on each electrical contact's node so the balance closes. |
| T6 | Thermal contacts | Device data on boundary patches, separate from electrical contacts: isothermal (own temperature) or thermal resistance R_th [K cm²/W, per boundary area in every dimension] to its own ambient temperature. Every other boundary is adiabatic. Contact temperatures are part of a bias point and can be swept. |
| T7 | Steady well-posedness | Steady, DC sweeps, `trace_bias`, AC operating points and quasi-static solves reject (`invalid_input`) a connected thermal domain with no isothermal or R_th thermal contact. Transient does **not** reject a fully adiabatic domain: with ρc ∂T/∂t and a defined initial T it is well posed, and the temperature evolves without a sink. |
| T8 | Traps and BTBT | Interface traps or band-to-band tunnelling combined with electrothermal return a typed `invalid_input` error. Not frozen at T_amb (that would be a silent physical approximation). Their electrothermal coupling is deferred to a follow-up unit. |
| T9 | Temperature paths with ∂/∂T | n_i via Varshni E_g(T) and Nc/Nv ∝ T^1.5; Caughey–Thomas mu_max; Boltzmann and Fermi–Dirac statistics; Slotboom BGN; the impact-ionization factor γ(T); ohmic contact densities and built-in potential at the contact T; gate φ_m(T) and the lumped-oxide gate rows; heterojunction band shifts and the thermionic emission velocity. SRH lifetimes, Auger coefficients and v_sat stay temperature-independent (legacy). |
| T10 | Temperature range | Newton steps in T are clipped (legacy 50 K, scaled by T_ref); a converged T outside the range `check_temperature` accepts is a typed error; the convergence criterion gains an absolute bound in T. |
| T11 | Thermal material data | κ(T) = κ300 (T/300)^-exponent and a constant ρc per material, every value cited (legacy: Si 1.48 W/cm K with 1.33, 4H-SiC 3.7). SiO₂ κ constant 0.014 W/cm K. A material with no thermal data used while electrothermal is on is an error; no values are invented. |
| T12 | Solve paths | All: steady bias and sweeps, `trace_bias` (thermal runaway folds), transient (ρc storage, BE/BDF2, LTE including T), AC (thermal storage in C). Equilibrium unchanged (T uniform at ambient). `Equations::equilibrium_poisson` with electrothermal is an error. |
| T13 | Off path | With electrothermal off every result is bit-identical to `main`. With it on at uniform T = T_ref, agreement to round-off (not bit identity). |
| T14 | Excluded | Hydrodynamic carrier temperatures, Kapitza interface resistance, radiation, temperature-dependent ρc, unstructured meshes (Unit 16). |

Gates (owner-approved): material laws and every ∂/∂T against FD, F_0/F_1 against 40-digit values; heat equation alone
(Kirchhoff-transform exact solution, Robin and adiabatic analytic, legacy G-PARABOLA and G-BC); FD Jacobian of the
four-block system 1D/2D/3D with the model set on; consistency (off bit-identical, κ×1e8 against the isothermal solve
within the ΔT bound, uniform 350 K against the isothermal device at 350 K); energy balance (steady exact, transient with
the storage terms); per-edge Joule part ≥ 0; Seebeck V_oc = ∫P dT read at the I = 0 crossing of a sweep (Boltzmann and
Fermi–Dirac); Peltier cooling at a current-carrying contact and a junction against analytic ΔT; zero currents and heat
at uniform-T equilibrium; self-heated resistor against a high-resolution C++ reference ODE (legacy G-ROLLOFF direction);
diode with R_th against the lumped model and the runaway fold located by `trace_bias` against the lumped criterion;
step-heated rod transient at orders 1 and 2; AC low-frequency limit equal to dI/dV of the electrothermal sweep, the
thermal pole on a fixture where the internal diffusion time is short, high-frequency limit equal to the isothermal
admittance; 2D/3D extrusion and MOSFET self-heating; run identity carries the thermal data; mutation checks and cost.

- **Owner decision:** [x] accept (2026-10-09)  [ ] change  [ ] defer

### Unit 23 implementation choices beyond T1–T14 (U1–U8, ruled on by the owner on 2026-10-09)

Each is recorded in `ARCHITECTURE.md` 6.2 "As built (Unit 23)". Rulings and what was done:

| # | Choice | Owner ruling | Implemented |
|---|---|---|---|
| U1 | The fourth unknown is the scaled rise τ = (T − T₀)/T₀, not θ = T/T₀. | **Accept**; correct the rationale; verify the derivatives and energy residuals. | Rationale corrected: a double near 1 is spaced 2.2e-16, so θ resolves T to 6.7e-14 K at 300 K; the sink heat is the flux κΔT/h through the edge next to the sink, resolved with θ to κ·6.7e-14 K/h (4e-8 W/cm² for h = 25 nm, not "across a micron"). Measured on a diode at 0.3 V: imbalance 9.0e-9 W/cm² with θ, 2.4e-12 with τ (the flows' rounding floor, U8; earlier wrongly attributed to the Newton tolerance). Tests: the heat rows' exact energy identity Σ K₀ f_T = −Σ V I at arbitrary states (to 1e-13); the FD Jacobian at rises of 1e-7 T₀ (1.5e-10). |
| U2 | The flux freezes T at the edge's mean θ_e in its coefficients (Scharfetter–Gummel in w = θ^(1+r) n). | **Accept**; add nonuniform-T convergence, equilibrium and Jacobian tests. | Tests: an open-circuit bar between 300 and 400 K carries no current on any edge (within its resolution) and its quasi-Fermi rise converges at order 2.00 by self-refinement (21 to 321 nodes), within 1e-5 of ∫P dT; the FD Jacobian across a factor-2.5 temperature jump on one edge, Boltzmann and Fermi–Dirac (4.4e-9). |
| U3 | Thermionic emission across a temperature step. | **Defer**; no claim without analytical justification and benchmarks; reject explicitly. | Electrothermal with a thermionic-emission interface is `invalid_input` ("DECISIONS.md U3, deferred"); the non-isothermal emission code is removed. Test: the refusal, and the same device accepted with the model off. |
| U4 | Kirchhoff transform in one material; the harmonic mean of the ends' κ at their own T between materials. | **Accept with gate**; add multilayer analytical and mesh-refinement tests. | Test: two layers in series (κ_A = 1.48 (T/300)^−1.33, κ_B = 0.3 (T/300)^0.5) between 300 and 500 K against the exact Kirchhoff solution: errors 4.5e-2 to 7.5e-4 K from 10 to 160 cells, orders 1.91 to 1.99. |
| U5 | Transient start. | **Change**: initialize per connected thermal domain; tell regions without a local sink from sink-free domains. | Each connected thermal domain with a thermal contact (a region with none of its own included, if joined to one) starts from its steady state; only domains with no thermal contact at all are held at T₀; a given initial temperature is held everywhere. Tests: a resistor with a sinkless oxide layer starts from the steady state everywhere; of two unconnected resistors, the one with a sink starts from its own steady state (exactly) and the other at exactly 300 K. |
| U6 | Incomplete ionization coupled with T. | **Change**: reject with `invalid_input` until the bound-state energy accounting is derived and validated. | Refused ("DECISIONS.md U6"); its temperature-coupled code (ionized doping, bound-carrier storage, contact terms) is removed. Test: the refusal, and the model off accepted. |
| U7 | `trace_bias` evaluates the bordered rows at the swept bias, dF/dλ at the state. | **Accept**; verify the full bordered Jacobian and the fold regression. | Test: the bias column against central differences of the residual in each contact's bias (3e-11), every nonzero on `bias_rows`, the R_th contact's heat row nonzero; with the FD-gated J this is the bordered matrix. The fold regression (4.6e-9 V from the lumped fold) is unchanged. |
| U8 | The steady energy-balance gate. | **Change**: separate absolute and relative gates; justify and test the threshold against the solver tolerances. | Two gates at every bias, at Newton tolerances 1e-8 and 1e-12: absolute \|heat − ΣVI\| ≤ E_max R_I (the two sides are the same energy flows summed two ways; their currents are resolved to R_I, `terminal_current_resolution`, E_max the largest carrier or contact energy in V), measured 1e-4 to 5e-3 of it; relative 1e-9 where E_max R_I < 1e-9 ΣVI (resolved points), measured ≤ 3.4e-12. The imbalance at 0.3–0.5 V (2.4e-12 W/cm²) does not depend on the tolerance; the tolerance shows at 0.7 V (1.4e-10 at 1e-8, 3.2e-11 at 1e-12). |

- **Owner decision:** U1, U2, U7 accept; U3 defer; U4 accept with gate; U5, U6, U8 change (2026-10-09).

---

## Commercial-readiness programme (owner request, 2026-10-09)

The owner asked for the architectural work a commercially viable device simulator needs, in this priority: input/output
and checkpoint/restart; realistic geometry and unstructured 3D meshing; robust nonlinear solving; low-current accuracy;
generic field/equation abstractions; scalable linear solvers; results and provenance; visualization, GUI and batch
workflows; validation against published and measured data. Constraints restated: foundations and numerical reliability
before more physics; Windows-only C++23, Win32/D3D12/D2D, no Python, no Qt; existing physics, APIs and validated
behaviour preserved; each change isolated and regression-tested in Debug and Release; no merges without authorization.
`ARCHITECTURE.md` 11.3 tracks each item's status (implemented and verified, in progress, planned).

### Unit 28 — bias step cutting (merged to `main` as `06b708d`, PR #31, 2026-10-14)

The branch starts at `59de5d5`, the head of `physics/electrothermal`, which PR #30 merged into `main` as `7300b67`
(2026-10-09); `main`'s tree equals `59de5d5`'s, so the branch carries Unit 28 alone against `main`.

Choices made without a ruling, open to the owner's review:

| # | Choice | Why | Trade-off |
|---|---|---|---|
| S1 | Off by default (`StepControl::max_cuts = 0`). | Every existing result, run identity and test stays bit-identical. | Users must turn it on; a later batch/input layer should default it on. |
| S2 | Halve from the last converged state; grow ×2 after a success; give up below 2^-max_cuts. | The standard commercial ramp; deterministic. | No predictor (the sub-step starts from the previous state with the new contact values stamped), so steps are shorter than an extrapolating ramp would need. |
| S3 | Only `non_convergence`, `singular_system`, `inaccurate_solve` are retried. | Cancellation and resource exhaustion are not cured by a shorter step. | — |
| S4 | Sub-step states are not sweep points; the point records `bias_substeps` and its last Newton run. | The sweep's points stay the requested ones, so extraction and comparisons are unchanged. | Intermediate states are lost (a checkpoint layer may want them). |
| S5 | Progress iterations count across a point's attempts. | Keeps the Unit 10 rule: events strictly increasing. | An event's iteration number is no longer the Newton iteration within one run. |
| S6 | The first point of a sweep from a given state is cut only if `initial_bias_V` says that state's bias. | The solver cannot infer the bias of an arbitrary state. | One more option to fill; a checkpoint format should store the bias with the state. |
| S7 | Step control settings and `initial_bias_V` enter the run record and identity only when on. | Runs without it keep their identity (as Units 19, 20, 23 did). | — |
| S8 | With the control on, every attempt (the full step included) is limited to `StepControl::max_iterations` (default 25). | Measured: without it each failed attempt ran 100 iterations and the MOSFET benchmark did not finish in 20 minutes; with it, 12 minutes. | A point that would converge slowly from the full step (more than 25 iterations) is reached by cutting instead: the same steady state, a different path. |

Finding (ARCHITECTURE.md 11.3): cutting rescues the MOSFET field-mobility drain sweep at the default damping, but
expensively (7 sub-steps per 0.1 V point). Unit 30 then measured a Newton line search against it (below).

### Unit 29 — the io layer (merged to `main` as `f69a7a8`, PR #33, 2026-10-14)

**Owner decision (2026-10-09):** "File I/O layer: Approved." D5's layer list gains `io` (depends on base, mesh,
results; solve does not use it). Choices made within that approval, open to review:

| # | Choice | Why | Trade-off |
|---|---|---|---|
| IO1 | One binary container (magic, version, kind, length, payload, FNV-1a 64) for every content. | One reader, one set of checks; exact doubles; fast. | Not human-readable; a text input deck is a separate, later unit. |
| IO2 | Keyed records: unknown keys skipped, absent keys defaulted. | New result fields need no version bump; old files keep loading. | Removing or retyping a field needs a format version change. |
| IO3 | FNV-1a 64 as the checksum. | Already the project's digest (run identity); detects corruption. | Not cryptographic: it detects accidents, not tampering. |
| IO4 | Atomic write by temporary file and rename. | A crash or full disk never leaves a half-written file under the target name. | Needs write access to the target's directory; no file locking between processes. |
| IO5 | A checkpoint is one converged point's state plus its bias, the mesh digest and the run identity. | Minimal, and enough to restart a sweep (with Unit 28's `initial_bias_V`). | Does not carry the device or options; the caller must restart with the same ones (the run identity lets it check). |
| IO6 | A checkpoint is refused on a mesh with another digest. | A state on another mesh is meaningless. | A remeshed device cannot restart from it (state transfer is Unit 17). |
| IO7 | Device descriptions and an input deck are not in this unit. | Each is a large surface (materials, contacts, interfaces, traps, thermal). | Restart needs the device rebuilt in code until then. |

---

### Unit 30 — Newton line search (merged to `main` as `f3e2811`, PR #32, 2026-10-14)

Owner request (2026-10-14): implement a Newton line search on a separate branch and benchmark the MOSFET case against
step cutting, convergence reliability and runtime, before deciding whether Unit 28 needs changes.

| # | Choice | Why | Trade-off |
|---|---|---|---|
| L1 | Backtrack along the damped step x + λ(x₁ − x), not along the raw correction. | ψ is clipped and densities bounded, so scaling the raw correction is undone by the clipping; convex combinations of x and x₁ keep every bound. | The direction is the damped one, not Newton's, so descent is not guaranteed. |
| L2 | Two tests: `residual` (Armijo on Σ F²) and `natural` (Deuflhard), natural the default when enabled. | Measured: residual 0/10 on the MOSFET, natural 10/10. | One more back-substitution per trial for natural (no factorization). |
| L3 | The natural test in the convergence test's measure (`update_size`). | It is the norm Newton already trusts; two alternatives measured worse (11.3). | That measure is erratic in reverse bias, so the search fails −5 V jumps plain Newton solves. |
| L4 | No trial passing: take the best, never λ = 0; convergence on the full correction. | The legacy line search reported λ = 0 steps as converged (6.2, Unit 21). | A hopeless iteration still moves (uphill if it must) and ends in non_convergence. |
| L5 | Off by default; settings in the run record only when on. | Every existing result and identity unchanged; the measured regressions make it unsafe as a default. | Users must opt in per problem. |
| L6 | An accepted λ = 1 trial's F and J are reused. | A search that always accepts λ = 1 costs nothing extra and is bit-identical to none. | — |

**Recommendation on Unit 28:** no change. **Owner decision (2026-10-14):** merged as is (PR #31, `06b708d`), S1–S8
accepted as they stand. Step cutting is the only strategy that converged in
every case measured (MOSFET and all six diode jumps); the line search is faster where it works (MOSFET 47 s against
540 s) but fails cases plain Newton solves. One tuning question remains open: with both on, S8's per-attempt limit of
25 makes step cutting cut points the line search would finish in 27–43 iterations (88 s against 47 s). Raising it
would speed that combination but slow step cutting alone (each failed attempt costs more), so it is left at 25 until
more cases are measured.

---

## Closed units: summaries (full records in `docs/decisions/`)

Each unit's complete record — the design rows, criteria, owner rulings in their own words, measurements and
as-built notes — was moved to its own file under [`docs/decisions/`](docs/decisions/) without changing a word (2026-11-14). Below
are the main points of each: scope, its decision IDs, outcome, the owner's rulings and accepted limitations, and where it
landed. Rule for later units: a unit's record stays here in full while it is open, then moves to `docs/decisions/` and
leaves a summary here. A citation elsewhere (code comments, ARCHITECTURE.md) of a row ID ("DECISIONS.md PV3", "DECISIONS
11.2", "DECISIONS Unit 55 15.13") is found in that unit's record file, unchanged.

### Meshing and geometry

- **Unit 16 — unstructured 2D and 3D meshing** ([record](docs/decisions/unit-16.md); `main` `a09cf43`, PR #35).
  - **Owner rulings R1–R4 (2026-10-14):** 2D and 3D in one unit; an in-tree generator with no new dependency;
    interfaces as node regions with straddled node pairs; analytic doping profiles.
  - **Constructions M1–M16:** exact predicates; a simplex import with checked duals (P0–P3, P5); a 2D generator
    (polygons, straddled interfaces, structured boundary layers); a 3D generator (convex polyhedra, contact patches,
    planar interfaces); P4-certified devices; meshes on disk.
  - **Evidence labels:** proven, cited, heuristic (capped), limitation.
  - **U14 (MOSFET gate):** unmet with isotropic triangles; met with a structured layer under the gate.
  - Every 2D and 3D gate passed.
- **Unit 49 — boundary fill between a structured layer's columns** ([record](docs/decisions/unit-49.md); `main` `5bb47c3`, PR #69).
  - Fixed a Unit 16 `generate_mesh` bug in the fill between layer columns.
  - CI's Release job failed at vcpkg checkout, so the merge was unverified on CI.
- **Unit 53 — FinFET geometry and meshing** ([record](docs/decisions/unit-53.md); `main` `a2a57a0`, PR #74).
  - **FD1–FD15, FC1–FC14:** non-convex regions (piecewise linear complexes, `extrude3d` sections); boundary recovery;
    junction lines and bends; patches cut by interfaces; the reduced and shaped FinFET decks (about 3e5 nodes,
    P4/P6-certified).
  - FC10 met: the shaped deck converges at all 17 points.
  - **Accepted as limitations (owner, 2026-10-31):** FC8 (flat gate, an error floor at the fin's top corners); FC9
    (13.9% low in on-current against the tensor reference, not diagnosed); FC14 (Debug +146 s over 120 s).
- **Unit 54 — 3D structured boundary layers in swept meshes** ([record](docs/decisions/unit-54.md); `main` `95c2254`, PR #75).
  - **UD1–UD13, UC1–UC13:** extrusions with 2D cross-section layers.
  - UC9 met: the shaped deck converges with 71,188 nodes, in 892 s.
  - **Accepted as limitations (owner, 2026-10-31):** UC7; UC8 (8,300 nodes, about as accurate as the coarse tensor
    mesh but outside the bound at 5 of 17 points).

### Solver robustness and numerics

- **Unit 36 — trace robustness at the current-resolution floor** ([record](docs/decisions/unit-36.md); `main` `07cd53a`, PR #46).
  - **R1–R3, TR1–TR7:** the arc row's products in compensated arithmetic; a corrector-distance test.
  - The zener band-to-band trace on the defaults now completes, and the open-base snap-back takes 59 points (was about 5000).
  - **Limitation:** the arc metric's current term still comes from the contact current.
- **Unit 41 — last-bit reproducibility on hybrid CPUs** ([record](docs/decisions/unit-41.md); `main` `67e8fe5`, PR #51).
  - Eigen's cache sizes are pinned to the P-core values, so results do not depend on core type.
  - Fix chosen by the owner (2026-10-17); the P-core results were kept.
- **Unit 44 — iterative linear solver** ([record](docs/decisions/unit-44.md); `main` `b518c99`, PR #54).
  - **LS1–LS6, LQ1–LQ4, LG1–LG8:** GMRES with ILUT, opt-in, real and complex, with a reported SparseLU fallback.
  - In 3D it is 9–72× faster and uses 8–26× less memory than SparseLU.
- **Unit 50 — Newton at the rounding floor** ([record](docs/decisions/unit-50.md); `main` `0c14841`, PR #71).
  - The residual line search takes the full step once the residual is at rounding (this unstopped the FinFET with MLDA).
  - Adds an opt-in `minority_floor` for uncontacted n+ regions under an oxide.
- **Parallel test suites** ([record](docs/decisions/parallel-test-suites.md); `main` `2fe1879`, PR #62).
  - Owner, 2026-10-23: "the test is taking too long", then "do the test-speed".
  - The suites run in parallel; the Release suite went from about 1,013 s to about 618 s.
  - The WG1 bound is 250 ms on CI.

### Physics models (each opt-in and bit-identical when off unless stated)

- **Unit 34 — Masetti hole mobility** ([record](docs/decisions/unit-34.md); `main` `e14a53f`, PR #42).
  - Criteria AC1–AC9; off by default.
  - The largest boron resistivity gap falls from 48% to 11%.
  - **Accepted as limitations (owner, 2026-10-15):** AC3 and AC5 at 1e18, as documented limitations, not as passed criteria.
- **Unit 35 — published silicon defaults** ([record](docs/decisions/unit-35.md); `main` `2d1fdc7`, PR #44).
  - A default change the owner authorized (S1–S5, DC1–DC10).
  - The band set is now Couderc et al. (2014): n_i(300 K) = 9.673e9.
  - The hole mobility is `automatic` (Masetti on silicon).
  - The legacy set and Caughey–Thomas stay selectable and reproduce old results byte for byte.
- **Unit 37 — density-dependent dopant ionization** ([record](docs/decisions/unit-37.md); `main` `16a77b5`, PR #47).
  - **I1–I5, II1–II7:** Altermatt et al.'s levels as `ionization_levels = altermatt`.
  - The incomplete-ionization gaps at 1e19 fall from +138%/+241% to +3%/+4%.
  - **Limitation:** donors are taken as phosphorus and acceptors as boron.
- **Unit 38 — Masetti electron mobility** ([record](docs/decisions/unit-38.md); `main` `1b446d0`, PR #48).
  - **E1–E4, ME1–ME7:** `electron_mobility = masetti`.
  - The largest phosphorus gap falls from 8.8% to 6.2%.
  - **Limitations:** μ_max = 1414 rests on one reproduction; phosphorus only.
- **Unit 46 — Lombardi surface mobility** ([record](docs/decisions/unit-46.md); `main` `51c2da9`, PR #64).
  - **LD1–LD8, LS1–LS12:** E⊥ per edge, an exact Jacobian; works with electrothermal; refused with impact ionization (LD8).
  - On the FinFET deck the on-current falls by 12%.
  - **Accepted as limitations of the model as specified:** LS6, LS8 and the mesh dependence.
- **Unit 47 — MLDA quantum confinement** ([record](docs/decisions/unit-47.md); `main` `329f203`, PR #66).
  - **QD1–QD10, QC1–QC13:** the density of states times 1 − exp(−d²/λ²); no new unknowns; both statistics.
  - Refused with electrothermal, thermionic interfaces and non-silicon regions.
  - On the MOS-C the peak moves 0.8 nm and V_T rises 26 mV.
- **Unit 48 — Schottky contacts** ([record](docs/decisions/unit-48.md); `main` `baa8a76`, PR #67; the ruling in PR #68).
  - SD1–SD9, SC1–SC17.
  - Merged while CI was running.
  - **Accepted as limitations (owner, 2026-10-28: "accept it as limitations"):** five unmet criteria, as limitations
    of the criteria as written, not defects of the model.
- **Unit 51 — gate tunnelling** ([record](docs/decisions/unit-51.md); `main` `cb3f9de`, PR #72).
  - **GD1–GD10, GC1–GC17, GR1–GR12:** Tsu–Esaki with WKB, through a lumped oxide or straight meshed paths; self-consistent with an exact Jacobian.
  - MLDA with tunnelling is refused (GR12).
  - **Limitation (GR5):** a tunnelling MOS-C with no minority-carrier contact is near-singular.
- **Unit 52 — energy-balance transport** ([record](docs/decisions/unit-52.md); `main` `85ad670`, PR #73).
  - **ED1–ED14, EC1–EC19:** both carrier temperatures as Newton unknowns; Scharfetter–Gummel energy fluxes.
  - The peak velocity in the ballistic diode is 1.93 v_s.
  - EC13's energy half was ruled "accept as built" (2026-10-29).
  - **Limitation:** with no E_g sink, T is overestimated near breakdown.
  - Many model combinations are refused.

### Inputs, outputs and results

- **Unit 31 — device files, the input deck and `nitcad`** ([record](docs/decisions/unit-31.md); `main` `6f67510`, PR #37).
  - **DK1–DK11:** our own sectioned `key = value` format with an in-tree parser (no TOML or JSON dependency); devices on disk; `nitcad check|run`.
- **Unit 32 — the weighted terminal current** ([record](docs/decisions/unit-32.md); `main` `6e54c69`, PR #39).
  - **W1–W8:** opt-in conductance-weighted current extraction with an adjoint error estimate.
  - The MOSFET off state and sub-floor diode currents are resolved.
- **Unit 42 — results and provenance** ([record](docs/decisions/unit-42.md); `main` `145f23c`, PR #52).
  - **PV1–PV7, PD1–PD5, PG1–PG8:** build, revision, platform, start and wall time in every run record and file; `nitcad --version`.
  - Accepted by the owner 2026-10-17.
- **Unit 43 — transient, small-signal and trace analyses in the deck** ([record](docs/decisions/unit-43.md); `main` `913dac0`, PR #53).
  - AN1–AN8, AQ1–AQ5, AG1–AG8.
  - Accepted by the owner 2026-10-17.
- **Unit 45 — rates and the mesh digest in results** ([record](docs/decisions/unit-45.md); `main` `81622ca`, PR #56).
  - RS1–RS6, RQ1–RQ4, RG1–RG8.
  - Accepted by the owner 2026-10-21.

### Architecture

- **Unit 40 — generic fields and equations behind the assembler** ([record](docs/decisions/unit-40.md); `main` `629bed0`, PR #50).
  - **FE1–FE7, FD1–FD9, FG1–FG9, FR1–FR9:** field descriptors and equation terms behind `DriftDiffusion`, bit-identical.
  - Deferred: the insulator nodes' unused unknowns (FD5); concurrent assembly (FD6); the "derived" field kind (FD7).
  - Accepted by the owner 2026-10-17.

### Validation against published data

- **Unit 33 — validation against published data** ([record](docs/decisions/unit-33.md); `main` `a14ac96`, PR #41).
  - **V1–V6:** references are values printed in the source's text or tables (no digitized curves), run as decks in `benchmarks/`.
  - Gaps are pinned, not tuned: boron mobility up to +48%, incomplete ionization at heavy doping, n_i +10%.
- **Unit 39 — breakdown validation** ([record](docs/decisions/unit-39.md)).
  - **Deferred** by the owner (2026-10-16): Sze and Gibbons (1966) and the reviews reproducing it are paywalled or behind a login.
- **Unit 55 — theoretical calculation validation framework** ([record](docs/decisions/unit-55.md); **deferred by the owner 2026-11-12**).
  - **Design (§1–§9):**
    - V&V levels L0–L5 (TVM1: no claim is made at a level until the levels below it pass);
    - reference classes E/D/A/M/H/P/X and independence grades I0–I3 (I0, a reference that calls NiTCAD code, is not evidence);
    - a `verify` layer whose reference library depends on `base` only (double-double arithmetic);
    - consistency checks C1–C12, catalogue entries TVB-01–34, design choices TVD-1–23;
    - questions TVQ1–6, ruled 2026-10-31 with three conditions (TVD-16–18).
  - **Rule TVM2:** never loosen a criterion, never tune a model.
  - **Phases:**
    - TVR-0, an audit of the existing validation (§10);
    - TVR-1, reference numerics (§11);
    - Unit 56, the fixes of findings (i)–(iii) (§12–§13);
    - TVR-2, code verification (§14);
    - TVR-3, the model-implementation catalogue (§15);
    - TVR-4, device theory (§16).
    - TVR-1 to TVR-4 are merged to `main` (PRs #76–#80).
  - **Open:** TVR-5a stopped after its step 5 (`verify/interfaces-generation`, `c8e36b8`, pushed, no PR; its later
    sections go into this record when it resumes). Still to be ruled: the S2 WKB first-order finding.

### Visualization and the Windows application

- **Units 25 and 26 — plan; Unit 25a, visualization** ([record](docs/decisions/units-25-26.md); 25a on `main` `6d2135b`, PR #55).
  - VR1–VR9, VQ1–VQ14, VG1–VG8, WG1–WG6.
  - CPU scenes; Direct3D 12 / Direct2D drawing to PNG; `nitcad plot`.
  - **Stack:** Win32 + D3D12 + D2D only.
  - Accepted by the owner 2026-10-21.
- **Unit 25b — derived fields, isosurfaces, streamlines** ([record](docs/decisions/unit-25b.md); `main` `8256cf8`, PR #57).
  - DS1–DS8, BQ1–BQ3, BG1–BG9.
- **Unit 25c — arrows, clipping box, exploded view, filmstrips** ([record](docs/decisions/unit-25c.md); PR #58).
  - DC1–DC6, CQ1–CQ4, CG1–CG9.
  - Accepted by the owner 2026-10-22.
- **Unit 26 — the Windows application** ([record](docs/decisions/unit-26.md); `main` `e9b1d71`, PR #59).
  - **AQ1–AQ4, WG1–WG8:** `nitcad-app.exe`, with a deck editor, a worker thread, a run monitor, results and `--self-test`.
  - Accepted by the owner 2026-10-22.
- **Unit 26b — the analysis workspace** ([record](docs/decisions/unit-26b.md); `main` `51face0`, PR #60).
  - **XD1–XD7, XG1–XG8:** live checking, structure preview, cut lines, run comparison, mesh quality.
  - The WG1 gap is deferred for the owner to investigate.
- **Unit 26c — GUI catch-up** ([record](docs/decisions/unit-26c.md); `main` `1be857e`, PR #61).
  - **YD1–YD6, YG1–YG6:** the 25b/25c views in the window, region visibility, playback, clip sliders.
  - The WG1 spikes are still open; the debug-layer self-test gap bound is 1000 ms (ruling 2026-10-31).

---

## Open units (records kept here in full)

### Unit 17 — intelligent adaptive mesh refinement and state transfer (proposed 2026-10-31; approved with amendments 2026-10-31; on `mesh/adaptive-refinement`)

**Owner request:** "I want Intelligent Adaptive mesh refinement" (2026-10-31). ARCHITECTURE.md 6.12.4 and the unit
table's Unit 17 ("Adaptive mesh refinement and state transfer", target). The branch started from `main` `95c2254`; rebuilt on `main` `03e0bf7` (PR #80) when implementation was authorized
(Unit 54 merged). The standing conditions hold: no change to any existing mesh, device, result or run identity when
adaptation is off; criteria fixed before code; Release and Debug suites at every commit, `/W4 /WX`; hand mutants;
nothing pushed or merged without authorization.

**Owner rulings (2026-10-31):** "Approve the proposed AQ1–AQ5 and LQ1–LQ4 recommendations. Proceed with the Unit 17 /
Unit 17b split." Amendments, incorporated below before the proposal was frozen:

| # | Amendment (owner) | Where |
|---|---|---|
| OA1 | "Treat greedy benefit-per-cost allocation as a heuristic, not a globally optimal solution." | ID6 (no optimality claim), AC16 |
| OA2 | "Require independent numerical error verification, conventional AMR fallback and explicit accuracy certification." | AD11–AD13, ID1, AC17, AC20 |
| OA3 | "Add false-positive refinement/mesh overhead measurements alongside the 80% regional recall target." | AC14 |
| OA4 | "Estimate computational cost using solver timings, DOFs and Jacobian sparsity—not node count alone." | ID4, ID5, ID6, AC15 |
| OA5 | "Keep cross-run ML opt-in, with reproducible training data, model versioning and held-out device validation." | ID8, ID10, AC19 |
| OA6 | "Keep numerical accuracy independent of ML predictions." | The governing rule below; AD11–AD13; AC21 |

"Do not write implementation code yet." The design and criteria below are frozen for review; implementation starts
only on the owner's word.

**Implementation authorized (owner, 2026-11-12):** "Implement Unit 17 adaptive refinement." The design and criteria
below are implemented as frozen (no criterion changed); the slices follow in order, each a local commit. The three
proposal commits (`d363963`, `9d879be`, `a534bda`, on the old base `95c2254`) were folded into one commit on
`main` `03e0bf7`; the text is unchanged apart from this paragraph and the base named above.

**Owner definition (2026-10-31):** "Intelligent Adaptive Mesh Refinement (Intelligent AMR): adds a predictive or
learning-based decision layer that can identify where refinement is likely to be needed, estimate its benefit, and
optimize mesh allocation." The unit therefore has two layers: the **adaptive base** (items 1–5 below: indicators,
size fields, transfer, the loop; a conventional AMR) and the **intelligent decision layer** on top of it (ID1–ID9,
below the base design): prediction, benefit estimation, allocation, and learning. **The governing rule (ID1, OA6):**
the decision layer only decides *where to spend nodes*. Numerical accuracy is independent of every prediction and
every learned model: the error estimate (AD3's adjoint), the independent verification (AD11), the stopping rule (AD6)
and the accuracy certificate (AD13) are computed from solves alone and never read a prediction, a fitted parameter of
the decision layer or a `.nitlearn` model. Every result is a converged solve on a certified mesh. A poor prediction
costs time, never accuracy; repeated poor predictions hand the run to conventional AMR (AD12).

**The adaptive base (the decision layer's foundation; each part gated below).** Today the user chooses every size: h,
per-region h, boxes, s, layers, tensor lines. The adaptive loop chooses them from the solution:
1. **Goal-driven.** The user states what must be accurate (a terminal current, or V_T / g_m / swing read from the
   sweep) and a tolerance, not mesh sizes. The loop stops when the goal's change between passes, and its estimate,
   fall below the tolerance, or at a node or pass budget, which it names.
2. **Physics-aware indicators.** Where the discretization is known to lose accuracy: the potential step across an
   edge in units of V_T (Scharfetter–Gummel's accuracy), the change of log density along an edge (inversion and
   depletion edges), doping steps (junctions, before any solve), and where active models concentrate (impact
   ionization and BTBT generation, the field at a tunnelling insulator).
3. **Goal-weighted (adjoint).** An indicator only counts where the goal is sensitive to it: each node's indicator is
   weighted by the adjoint z of the goal, J^T z = ∂I/∂u (one transposed solve per goal per adaptation point). A
   junction far from the current path (a well, the substrate) is then not refined for a drain-current goal.
4. **Bias-aware.** The user lists the bias points to adapt at (or `auto`: the sweep's ends and the point of largest
   goal curvature from a first coarse sweep); marks from every point are united, so one mesh serves the whole sweep.
5. **Explains itself.** A report per pass: nodes, the goal and its estimate, which indicator marked how many nodes
   and where (regions), the stopping reason; the final mesh saved and viewable like any generated mesh.

**The questions that shape the design (owner rulings 2026-10-31: every recommendation approved):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| AQ1 | How the mesh changes | (a) **Regenerate** with a size field: the existing generators already size every element through one `size_at(point)` (`generate2d.cpp:345`, `generate3d.cpp:554`, `generate_complex.cpp:785`); the indicators give a nodal h on the old mesh, graded, and the generator meets it. Every certificate (P0–P5, P4, P6) comes from the existing checked path. (b) **Local insertion** into the existing triangulation (Ruppert/Bowyer–Watson points where marked). Cheaper per pass, but P4's interface bands and the layers' exemptions must be re-proven for incremental insertion. | **(a)**: no new geometry, certificates unchanged. Cost is a full mesh per pass (2D: seconds; swept: the 2D cost; general 3D: about 130–160 s per 1e5 nodes, FC12). |
| AQ2 | Mesh kinds in this unit | (a) **Tensor** (insert lines where marked, per axis; nonlocal BTBT stays valid), **2D generated**, **swept 3D** (adapt the cross-section's size field and the axis lines). (b) Also general 3D (`generate3d`, Unit 53's complex). | **(a)**. All three are the 2D generator plus line sets; general 3D re-mesh cost is high and goes to a later unit. Structured layers (U14, UD5) are kept as given: never refined inside; an indicator saturating inside a layer is reported. |
| AQ3 | Mesh per sweep or per bias point | (a) **One mesh per sweep**: adapt at the listed points, refine-only union of marks, then the whole sweep on the final mesh. (b) A mesh per point (state transferred point to point). | **(a)**: one mesh, one run identity, curves comparable point to point; sweeps, files, plots and the app unchanged. (b) changes every downstream layer. |
| AQ4 | Goal weighting (adjoint) | (a) **In this unit**: a transposed solve in `linalg` (SparseLU and PARDISO; GMRES-ILUT falls back to SparseLU for this solve, reported) and ∂I/∂u from the weighted current (Unit 32's W1, which is linear in the state's fluxes). (b) Feature indicators only; adjoint later. | **(a)**: it is what makes the loop goal-driven rather than "refine everything steep". Gated by AC7 against feature-only. |
| AQ5 | Coarsening and anisotropy | Coarsening below the initial mesh's sizes; anisotropic (metric) refinement. | **Both out of this unit.** The initial mesh is the coarsest; anisotropy is already supplied by layers and tensor lines, which the loop adapts along their own axes (tensor lines, swept axis lines). |
| LQ1 | What "learning" learns from | (a) **Within a run** (online): the decision layer fits its own models to the run's passes and bias points (local convergence orders, error constants, cost per node) and uses them for the next decision. (b) **Across runs** (persistent): a model trained on earlier runs' adaptation histories (a versioned, checksummed `.nitlearn` file, its digest in the run identity) predicts the first mesh of a new run. (c) **Shipped pre-trained** model: needs a training corpus that does not exist yet. | **(a) and (b) in this unit, (b) opt-in**; (c) not now. (a) needs no stored state and is fully reproducible from the run. (b) is where "likely to be needed before solving" pays off most (fewer passes from the first mesh). |
| LQ2 | Model class for the cross-run predictor | (a) **Gradient-boosted regression trees**, in-tree C++ (no ML library; no Python, owner rule), deterministic training. (b) Ridge regression on engineered features. (c) A small neural network (MLP), in-tree C++. | **(a)**: handles mixed features (material, region kind, log doping, distances, bias) and their thresholds without scaling, deterministic, small (hundreds of trees × depth 4), inspectable (feature importances in the report). (b) is the baseline it must beat (ID9). (c) needs far more data than a user's run history gives. |
| LQ3 | The allocation objective | (a) **Tolerance mode**: reach the goal tolerance with the fewest nodes. (b) **Budget mode**: the smallest estimated error for a node (or wall-time) budget. | **Both**, the user picks; tolerance mode the default. |
| LQ4 | One unit or two | (a) One unit, all of it (about 9 slices, XXL). (b) **Unit 17: base + online decision layer (prediction from physics and from the sweep, benefit, allocation; LQ1 a); Unit 17b: the cross-run learned predictor (LQ1 b, LQ2).** | **(b)**: 17b needs 17's adaptation histories as its training data, and 17's gates give the baseline 17b must beat. The design of both is fixed now. |

**Design (as approved, with amendments OA1–OA6):**

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| AD1 | **The loop** (`solve::adapt`, new): solve (the existing analysis) → indicators → mark → size field → regenerate → transfer → solve; the analysis and the device builder are passed in, so physics, solver and analysis code do not change. | 6.12.4's target loop. | Device-neutral (6.12.8). | One more layer of orchestration in `solve`. |
| AD2 | **Indicators** per edge, each normalized to a dimensionless "excess" e ≥ 0 (0 = resolved): Δψ / (k V_T) with k = 1 default; \|Δ ln n\|, \|Δ ln p\| / ln 10 per decade; \|Δ ln \|N_D − N_A\|\| for doping; generation-rate integral fraction for II/BTBT (only when the model is on). A node's indicator is the largest over its edges. | SG accuracy degrades with Δψ ≫ V_T (standard); the rest heuristic, measured by AC5–AC7. | Cheap, device-neutral, computed from fields the solver already has. | Thresholds are parameters with defaults; the defaults are fixed in slice 2 from AC4 measurements and recorded. |
| AD3 | **Goal weighting:** node weight ω_i = \|z_i\| scaled by the node's share of the goal (normalized so Σ ω_i = 1 over semiconductor nodes); marked if e_i ω_i ≥ θ max(e ω) (Dörfler-style bulk marking, θ = 0.3 default). Feature-only mode: ω ≡ 1. | Dörfler marking: standard in adaptive FEM; the adjoint weight: dual-weighted residual idea, here as a weight, not a bound. | Concentrates nodes where the goal is sensitive. | An estimate, not a bound; AC7 measures it. |
| AD4 | **Size field:** h_new = h_old / 2 at marked nodes (per pass), h_old elsewhere; then graded so h changes by at most a factor g = 1.5 per edge (fixed-point sweep), and never below `min_size` (default s / 4 where interfaces exist, else h_initial / 64). The generator queries it by exact point location on the old mesh (barycentric, linear). Tensor: a line inserted at each marked cell's midpoint along each axis where that axis' indicator marks it. | Grading: Unit 49's sliver finding. | Fed into the generators' existing `size_at`. | A new option field per generator (`size_field`); absent = bit-identical. |
| AD5 | **State transfer:** ψ, the quasi-Fermi potentials (not densities: densities interpolate badly across decades), T and carrier temperatures when present, linearly interpolated at the new nodes by exact location; interface and contact quantities recomputed. Used as Newton's initial state only: the result on the new mesh is a converged solve, so every conservation gate applies unchanged. Also exposed as `transfer_state` so a checkpoint can restart on a remeshed device (IO6's gap), still refused by digest unless transferred explicitly. | — | Fewer Newton iterations per pass; checkpoints usable after remeshing. | Interpolation across a straddled interface uses the node's own side. |
| AD6 | **Stopping:** the loop proposes to stop when \|G_k − G_{k−1}\| ≤ tol·\|G_k\| for every goal and the adjoint estimate ≤ tol; the stop is then subject to AD11's verification. Otherwise it ends at `max_passes` (default 6), `max_nodes` or a DOF/time budget, uncertified (AD13). The reason recorded. | — | The user's stated accuracy, not a mesh size. | Two consecutive small changes can be stagnation; the estimate and AD11 guard it. |
| AD7 | **Deck:** `[adapt]` section: `goal = current <contact>` (repeatable) or `vt`, `gm`, `swing`; `tolerance`; `at = <bias points>` or `auto`; `indicators = auto` or a list; `weighting = adjoint` or `none`; `policy = intelligent` or `basic`; `mode = tolerance` or `budget` with `budget_dofs` / `budget_seconds`; `verify = richardson` or `off`; `certify = required` or `report`; `max_passes`, `max_nodes`, `min_size`, `grading`, `theta`; (17b) `predictor = <file>`. Absent: no adaptation, bit-identical. | DK conventions. | The batch path reaches it without code. | — |
| AD8 | **Results, io, run record, app:** the adaptation report (per pass: nodes, goals, estimate, marks by indicator and region, time, stopping reason) in results and the run record; the final mesh saved as a generated mesh with its digest; the run identity covers the `[adapt]` settings and the final mesh digest; the app shows the passes' meshes and the report. | Unit 42 provenance rules. | Reproducible and explainable. | — |
| AD9 | **Determinism:** marking ties broken by node index; the size field and regeneration are deterministic; two runs bit-identical on one machine and thread count (6.8). | — | — | — |
| AD10 | **Existing paths:** with no `[adapt]` section and no `size_field`, every generator, device, result, file and run identity is bit-identical. | — | — | — |
| AD11 | **Independent numerical error verification (OA2):** before a result is certified, the final mesh M is refined uniformly once (every size field value halved, every tensor and axis line interval split; layers' lines halved) to M′, the analysis is solved on M′ from the transferred state, and the verified error of each goal is the Richardson estimate \|G_{M′} − G_M\| / (2^p − 1) with p the *observed* order where a third mesh exists from the loop's passes (clamped to [0.5, 3]), else p = 1 (conservative for the box method at abrupt features). This uses neither the indicators, nor the adjoint, nor any prediction: a second, independent check of AD3's estimate. The solve on M′ is reported but is not the result (the result stays on M, one mesh per sweep). `verify = richardson` is the default; `verify = off` skips it and the result is marked unverified, never certified. | Richardson extrapolation (standard); its independence from the adjoint estimate is by construction. | An adjoint estimate can be wrong where the linearization is poor (strong injection, near breakdown); a refined solve is not linearized. | Cost: one solve on about 4× (2D) or 4–8× (swept) the DOFs, measured (AC11). |
| AD12 | **Conventional AMR fallback (OA2):** the run switches from `policy = intelligent` to `policy = basic` (AD3's marking, no predictions) for every remaining pass when (a) two consecutive proposals are misses (ID1), (b) a pass lowers neither the estimate nor the goal change, or (c) AD11 rejects a proposed stop (verified error > tol): the next pass then marks by AD3 using the verified error's distribution. The switch, its trigger and its pass are recorded. | — | The intelligent layer can only speed up a conventional AMR, never replace its guarantees. | A run that falls back is no faster than basic, but no less accurate. |
| AD13 | **Explicit accuracy certification (OA2):** every adaptive result carries, per goal: the value, the adjoint estimate (AD3), the verified error (AD11) and the order used, the tolerance, the policy that produced the mesh and whether a fallback occurred, and a status: **certified** (estimate ≤ tol *and* verified error ≤ tol), **not certified** (either exceeds tol, with the margin) or **unverified** (`verify = off`, or a budget stop). The status, numbers and method are in results, the run record, the CSV table and the app; `nitcad run` exits with a distinct non-zero code when a requested `certify = required` run ends uncertified. | — | A user sees, per result, whether the stated accuracy was reached and how it was checked. | A certificate is evidence from two estimates, not a proof; the text says so. |

**The intelligent decision layer (`solve::adapt_policy`; replaces AD3's fixed-θ marking when on, `[adapt] policy =
intelligent`, the default once AC14–AC17 pass; `policy = basic` keeps AD3):**

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| ID1 | **Propose, then verify.** The layer proposes the next size field; the base loop builds, transfers, solves and estimates. Stopping (AD6), verification (AD11) and certification (AD13) never use a prediction; only solves. A proposal that does not lower the estimate is recorded as a miss and AD3's marking is used for that pass; two consecutive misses hand the run to conventional AMR (AD12). | — | Predictions can be wrong; results cannot. | A miss costs one pass. |
| ID2 | **A priori prediction (before any solve):** a predicted size field from the device alone. Junctions: depletion width W = sqrt(2 ε (V_bi + \|V\|) / (q N_eff)) at the largest bias of the analysis, from the sampled doping (M15); the Debye length L_D = sqrt(ε V_T / (q N)) sets h ≤ L_D / 2 across each junction's predicted span. Insulator interfaces under gates: the inversion-layer scale from the gate's largest bias and oxide capacitance (classical; MLDA's λ when on). Contacts: none. Active models: II and BTBT where the predicted junction field exceeds the model's onset. | Textbook depletion and Debye scales; heuristic as a mesh rule. | Pass 0 starts near the answer instead of from uniform h: fewer passes. | A junction whose bias the analysis does not reach is over-predicted; ID5 learns it down. |
| ID3 | **Prediction across bias:** with one mesh per sweep (AQ3), the layer solves the sweep coarsely once, tracks each node's indicator over bias, and chooses the adaptation points where the indicator field changes most (the largest change in its weighted norm between neighbours), up to `max_points`; features that move with bias (pinch-off, depletion edges) are covered by the union of their predicted positions. | — | "auto" adaptation points chosen from the physics instead of the sweep's ends. | One coarse sweep of cost. |
| ID4 | **Benefit estimate per candidate:** candidates are graded patches (a marked node and its neighbours to the grading radius), tensor lines, or swept axis lines. Predicted error reduction ΔE_c = Σ_{i∈c} ω_i η_i (1 − 2^{−p_c}), η the node's indicator excess (AD2) in the goal's units through the adjoint (AD3), p_c the local convergence order. **Predicted cost (OA4), not node count alone:** C_c = ΔT_pass(DOF + ΔDOF_c, nnz + Δnnz_c) − ΔT_pass(DOF, nnz), with ΔDOF_c = added nodes × the equation set's unknowns per node (3 drift-diffusion, 4 electrothermal, 5 energy balance; insulator nodes 1), Δnnz_c the Jacobian nonzeros the candidate adds (added edges × stride², counted from the candidate's predicted local mesh: mean degree from the current mesh's region, 2D and swept separately), and ΔT_pass the cost model of ID5. | The (1 − 2^{−p}) factor: halving h with local order p (standard a priori estimate). Sparse direct cost grows faster than DOFs (fill-in), so DOFs and nonzeros, not nodes, predict time. | Benefit in the goal's units, cost in seconds: comparable across candidates and across equation sets. | p and the cost model are learned (ID5); before they are, p = 1 (the box method's order at abrupt features, measured in Units 16 and 52) and the cost model's defaults from AC11's measurements. |
| ID5 | **Online learning (LQ1 a):** after each pass the realized error change per patch (the estimate before and after, through the adjoint) gives a fitted local order p per indicator kind and region (least squares in log2, clamped to [0.5, 3]). **The cost model (OA4)** is fitted from the passes' own measured timings: T_pass = N_newton × (a · nnz(J)^b + c · nnz(L+U) + d · DOF) + T_mesh(nodes) + T_transfer, with the factorization's fill nnz(L+U) as reported by the backend (SparseLU and PARDISO report it; GMRES-ILUT reports the preconditioner's), N_newton the measured Newton iterations per point, and T_mesh, T_transfer measured; b clamped to [1, 2]. Each fit's sample count, residual and predicted-vs-measured pass time go in the report. | — | The benefit and cost estimates become the device's and machine's own, not defaults. | Few passes, few samples: fits regularized to the defaults (ridge toward p = 1, b = 1.5). The decision layer only; AD11 and AD13 never read them (OA6). |
| ID6 | **Allocation (LQ3; a heuristic, OA1):** tolerance mode: take candidates in decreasing ΔE / C until the predicted error ≤ tol / 2 (margin for misprediction) or the pass' DOF growth cap (×2.5 per pass); budget mode: the same order until the DOF or time budget is spent. Ties by node index. **No optimality is claimed:** greedy by ratio is optimal only for a fractional knapsack with fixed, additive values; here the values are predictions, overlapping patches and grading couple the candidates, and refining one region changes the error elsewhere, so the allocation is a heuristic. | Fractional knapsack's greedy optimality is the motivation only, not a property of this problem. | A simple, explainable ordering that is measured against its baseline. | Can be beaten by a better heuristic; AC16 and ID9 measure it against AD3 and record the gap. |
| ID7 | **Explanation:** per pass the report lists the top candidates (location, region, indicator, ΔE predicted vs realized, cost), the fitted p and cost model, misses (ID1), and the predicted vs actual node count and error. | — | A user can see why the mesh is where it is. | Larger reports (capped at the top 20 per pass). |
| ID8 | **Cross-run predictor (Unit 17b if LQ4 b):** training samples are (node features → realized benefit) from saved adaptation histories (every run with `[adapt]` can save its history, opt-in). Features: material, region kind, log \|N_D − N_A\|, \|∇ log N\|, distance to the nearest junction, insulator interface and contact, predicted W and L_D (ID2), the analysis' bias range, the goal kind. Model: gradient-boosted trees (LQ2), fixed hyperparameters, deterministic training (sorted samples, no randomness); `[adapt] predictor = <file>` replaces ID2's a priori field with the model's, ID1 and AD12 still guarding. **Opt-in only (OA5):** off unless `predictor` is given; never on by default, never shipped. **Reproducible training data (OA5):** an adaptation history file (versioned, checksummed, atomic, as the `io` layer's files) holds the deck text, the device and mesh digests of every pass, the provenance (Unit 42), the per-node features and realized benefits, and the feature-schema version; `nitcad learn` builds a training set from a list of history files into a manifest (each file's digest, sorted), so the same manifest always gives the same samples. | — | Learns the user's device families (a MOSFET process with varying L, a diode series). | Trained on the user's own runs only; extrapolation outside them falls back to ID2 (features out of the training range are flagged). |
| ID9 | **The layer must beat its baseline:** every intelligence gate (AC14–AC17) compares against `policy = basic` (AD3), and 17b's against 17's online layer and against ridge regression; a component that does not beat its baseline is not made the default. | — | "Intelligent" is measured, not asserted. | — |
| ID10 | **Model versioning and held-out validation (17b, OA5):** a `.nitlearn` file records its format version, the feature-schema version, the model class and hyperparameters, the training manifest's digest, the build that trained it (Unit 42 provenance), and its held-out metrics. **Held-out devices, not held-out nodes:** the training command splits by device (history file), never by node, so no device is in both sets; the metrics (rank correlation of predicted vs realized benefit, recall/precision of predicted refinement, passes saved) are measured on the held-out devices and stored in the file. A model whose schema or format version differs from the build's is refused (`invalid_input`, both versions named); the model's digest and version enter the run identity and the run record. | Device-level splitting: nodes of one device are strongly correlated, so a node split overstates accuracy (standard leakage argument). | Any prediction a run used can be traced to the exact model, data and code that made it. | Retraining after a feature change is required; old models are refused, not silently reused. |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| AC1 | Nothing changes: every existing test and pin passes; every existing deck's run identity, mesh digests and device files byte-identical. |
| AC2 | Size field honoured: on the 2D generator, swept and tensor paths, every element's size ≤ the field at its centroid (the generator's own bound), grading ≤ g measured, certificates P0–P5, P4, P6 as for any generated mesh; a constant field equal to h gives the h mesh's digest. |
| AC3 | Transfer: a linear ψ and linear quasi-Fermi potentials transfer exactly (1e-12); a diode and a MOSFET state transferred to a refined mesh converge to the same solution as a solve from scratch (1e-9 relative in currents) in fewer Newton iterations (measured, reported). |
| AC4 | Indicators: on a 1D-in-2D abrupt diode the doping and Δψ indicators mark the junction's depletion edges only; on a MOS-C in inversion the Δ ln n indicator marks the inversion layer; each indicator's edge values equal hand-computed ones on a 4-node fixture. |
| AC5 | **Efficiency gate (2D):** the 2D MOSFET (Unit 16's deck shape) from a coarse mesh, goal = I_D at the sweep's points: the adaptive mesh reaches \|I_D − I_ref\| ≤ 1% at every point above the current resolution (I_ref from a 4× uniformly refined reference) with at most 1/3 of the nodes the uniform-refinement ladder needs for the same error. Risk stated before the run: if it fails, `[!shouldfail]` with its measurement for an owner ruling; never loosened. |
| AC6 | **Diode gate:** the legacy diode's forward and reverse currents reach the 2× reference within 0.1% with fewer nodes than the existing hand-made mesh; the loop stops by AD6's goal rule, not a budget. |
| AC7 | **Goal weighting earns its place:** a MOSFET with a deep well junction off the current path; adjoint-weighted vs feature-only to the same I_D error: the weighted mesh has ≥ 25% fewer nodes, and its marks avoid the well. Risk stated: if it fails, `[!shouldfail]`, owner ruling. The adjoint estimate is within a factor 3 of the actual error at the final pass. |
| AC8 | **FinFET (swept):** from a coarse swept deck, adapt for I_D; report nodes and FC9/UC8's rule against the 2× tensor reference at every point. Measured target: UC8's rule at all 17 points within UC8's 9,367-node budget; stated before the run, not a pass/fail gate (UC8 failed 5/17 by hand). |
| AC9 | Deck, io, run record, app: `nitcad check` validates `[adapt]` with line numbers; the report and final mesh round-trip; the run identity covers settings and mesh digest; `nitcad plot` and the app show the passes. |
| AC10 | Determinism: two runs bit-identical (AD9); refusals (unknown goal, contact, indicator; goal-weighting with a backend that cannot transpose and no fallback; adapt on general 3D in this unit) with code and message. |
| AC11 | Cost: per pass, DOFs, Jacobian nonzeros, factorization fill, Newton iterations, and regeneration + transfer + solve time reported for AC5, AC6, AC8; AD11's verification solve reported separately (time and DOFs); the transposed solve ≤ the cost of one Newton iteration's linear solve. |
| AC12 | Hand mutants: at least one per design row (AD1–AD13, ID1–ID7, ID9; 17b: ID8, ID10), each survivor closed by a test or argued equivalent. |
| AC13 | Release and Debug suites at every commit, `/W4 /WX`; Debug grows by ≤ 60 s; larger cases `[.benchmark]`. |
| AC14 | **Prediction (ID2, ID3), recall and false positives (OA3):** on the diode, the MOS-C and the 2D MOSFET, measured against the region the basic loop finally refines (nodes whose final size is below the initial h): **recall** (the a priori field covers ≥ 80% of it); **false-positive refinement** reported for every device: the precision (the fraction of the predicted-refined region that the basic loop also refined; its complement is the false-positive fraction), the false-positive DOFs and nonzeros, and their cost (the time of the final pass with and without them, the latter by re-solving on the basic loop's mesh), gated: false-positive DOFs ≤ 30% of the intelligent final mesh's DOFs; **mesh overhead**: the intelligent final mesh's DOFs ≤ 1.1 × basic's at the same verified accuracy. And the intelligent loop reaches AC5's and AC6's accuracy in at least one pass fewer than `policy = basic`, summed over the three. Risk stated: `[!shouldfail]` and owner ruling if not; never loosened. |
| AC15 | **Benefit and cost estimates (ID4, ID5, OA4):** over all candidates refined in AC5's run, the predicted and realized error reductions have Spearman rank correlation ≥ 0.6; the fitted orders on the abrupt diode lie within [0.8, 1.3] (first order expected); from the third pass on, the cost model predicts each pass' measured time within ±30% on AC5, AC6 and AC8, and does so better than a node-count-only model fitted to the same passes (mean absolute relative error, both reported). |
| AC16 | **Allocation (ID6, a heuristic, OA1):** tolerance mode: the intelligent loop's final DOFs ≤ 0.8 × `policy = basic`'s at the same verified accuracy on AC5's MOSFET and AC7's well MOSFET; budget mode: at a DOF budget of 50% of basic's final DOFs, the intelligent loop's verified error is ≤ basic's at that budget. Results stated as measured comparisons, with no optimality claim. Risk stated as AC14. |
| AC17 | **The guard and the fallback (ID1, AD12):** a deliberately wrong prediction (a hand-written field refining only the substrate) still reaches the tolerance and is certified, with the misses reported, the switch to `policy = basic` recorded at the pass AD12 names, and the result within the tolerance of the reference; total time ≤ 2 × basic's. Each AD12 trigger (a, b, c) is forced once in a fixture and its switch recorded. |
| AC18 | Determinism and explanation: two runs bit-identical, report and certificate included; the report's predicted vs realized numbers equal those recomputed from the saved passes. |
| AC19 (17b) | **Cross-run learning (OA5):** a family of 8 MOSFETs (L, oxide and doping varied, generated decks); held-out by device: trained on 6, tested on the 2 held out: the learned first mesh needs ≥ 1 pass fewer than ID2's a priori field and ends with ≤ the DOFs at the same verified accuracy; gradient-boosted trees beat ridge regression on held-out benefit prediction (rank correlation); the training set built twice from the same manifest is bit-identical, and so is the trained model; a model with another schema or format version is refused; the run identity changes with the model's digest; an out-of-range device is flagged and falls back to ID2; with `predictor` absent every run is bit-identical to Unit 17's. |
| AC20 | **Certification (AD11, AD13, OA2):** on AC5, AC6, AC7 and AC8 every adaptive result carries a certificate; against the references (AC5's 4× uniform, AC6's 2×, AC8's 2× tensor) no result marked **certified** has an actual error above its tolerance; a run made to stop early (a `max_nodes` budget) is marked **unverified** or **not certified**, never certified; with `certify = required` such a run exits with the distinct code; the verified error and the adjoint estimate both reported, and their ratio to the actual error recorded. Mutants on AD11 (skipped refinement, wrong p, wrong mesh compared) and AD13 (status logic) are killed. |
| AC21 | **Accuracy independent of predictions (OA6):** for AC5's and AC6's final meshes, the estimate, verified error and certificate computed with `policy = basic`, `policy = intelligent` and (17b) a `predictor` given are bit-identical when the final mesh is the same (a fixture forcing the same mesh), and a mutant making AD11 or AD13 read any decision-layer quantity is killed by this test. |

**Slices** (each a local commit with both suites green):
1. Size fields in the 2D generator, swept and tensor paths; exact point location on the old mesh; grading
   (AD4, AD10; AC1, AC2).
2. Indicators and marking, feature-only (AD2, AD3 with ω ≡ 1; AC4).
3. State transfer, the loop with stopping, independent verification and the certificate; the diode gate
   (AD1, AD5, AD6, AD9, AD11, AD13; AC3, AC6, AC20 part).
4. Transposed solve and goal weighting (AD3; AC7).
5. The decision layer: a priori and bias prediction, benefit, the DOF/sparsity/timing cost model, online fits,
   allocation, guard and conventional fallback, explanation (ID1–ID7, ID9, AD12; AC14–AC18, AC21).
6. Deck, results, io, run record, app, certificate display and exit code (AD7, AD8, AD13; AC9, AC10).
7. Gates and cost: AC5, AC8, AC11, AC20, mutants AC12, ARCHITECTURE 6.12.4 and 11.3.

Unit 17b (approved split, LQ4 b; a separate branch and unit after 17 merges): adaptation-history files and
`nitcad learn` manifests, features and schema versioning, deterministic gradient-boosted trees and the ridge
baseline, device-level held-out validation, versioned `.nitlearn` io, deck `predictor` (opt-in), AC19 and AC21's
predictor case (ID8, ID10); 3 slices.

**Size:** Unit 17 XL (7 slices; 1, 3, 4 and 5 the bulk); Unit 17b L.

**Implementation record (from 2026-11-12).**

*Slice 1 — size fields (AD4, AD10; AC1, AC2).* `mesh/size_field.hpp`: `SizeField` (a nodal size on a background
triangulation; exact location by `orient2d` through a bucket grid, the lowest-indexed containing triangle on a shared
edge, linear interpolation clamped to the corners' range; outside every triangle the value at the nearest point of the
nearest triangle), `local_sizes` (a node's longest incident edge, the h_old of AD4), `grade_sizes`, `refine_axis`.
`MeshOptions2D::size_field` (a `shared_ptr<const SizeField>`, null by default) enters the 2D generator only through
`size_at` (`min` with the field), so a field only refines; the swept generator passes its section's field through and
adds the field's digest to `generation_digest` only when one is given (AD10: no word added without it). Two
construction points found while writing the tests, fixed before any gate run (no criterion touched):
- *Grading follows the refinement only.* Grading every node of the old mesh (AD4's "fixed-point sweep" read
  literally) spreads the s-scale sizes of every interface pair and straddle over the whole device, refining where
  nothing was marked. `grade_sizes` therefore takes seed nodes (the nodes a pass refined) and propagates g^k h from
  them only (Dijkstra on log h); with every node a seed it is the full grading (h_j ≤ g h_i on every edge, tested).
  `refine_axis` likewise grades against intervals split in the call, keeping the given axis' own ratios.
- *Axis grading factor at least 2.* A midpoint split leaves a ratio of exactly 2 beside it, so with g = 1.5 every split
  cascades over the whole axis. On axes the factor is max(g, 2) (the 2D field keeps g = 1.5). Stated in the header.
Tests (`tests/mesh/size_field_test.cpp`, 6 cases): linear reproduction (1e-12) and the outside rule; refusals;
local sizes and grading (all seeds, seeds only); axis refinement (marked, second-level grading, protected straddles,
min_interval, refusals, determinism); AC2 on the 2D generator (a constant field equal to h gives the h mesh's
`graph_digest` and points; a left-half h/4 graded field: P4 certificate kept, more nodes, concentrated where the field
is small, every triangle not exempt within the generator's bound circumradius ≤ size/√3 against the field at its
centroid, deterministic); AC2 on the swept path (same digest of the graph with a non-refining field, generation digest
names the field, a refining field adds nodes in whole axis lines).

*Slice 2 — indicators and feature-only marking (AD2, AD3 with ω ≡ 1; AC4).* `solve/adapt_indicators.hpp`:
`edge_indicators` (raw values per edge and kind), `compute_indicators` (excess e = max(0, r/t − 1), node max, dominant
kind), `unite` (AQ3's union over adaptation points), `mark_nodes` (Dörfler: e w > 0 and e w ≥ θ max). Insulator edges
carry no indicator (the potential is linear between charges there, exact for the box method). The generation indicator
is a node's share of the total G·V times the number of semiconductor nodes, t = 10.

| # | Question (put and ruled under the owner's standing instruction "always choose recommended", 2026-11-11) | Options | Recommendation and ruling |
|---|---|---|---|
| AR1 | AD2's doping indicator \|Δ ln \|N_D − N_A\|\| is near zero across an abrupt junction between equal and opposite dopings (1e17 to −1e17: Δ ln = 0), so AC4's "the doping indicator marks the junction" cannot be met as written. | (a) As written. (b) **The signed log L(N) = sign(N) ln(max(\|N\|, n_i)/n_i), indicator \|ΔL\| / ln 10**: on an edge of one sign above n_i it is AD2's formula exactly; across a junction it counts the decades down to n_i and back up. (c) The equilibrium potential step \|Δ asinh(N / 2n_i)\|. | **(b)**, ruled as recommended: it keeps AD2's meaning and units (decades of net doping) where AD2's formula works and needs no solve; (c) changes the unit. No criterion changed. |

AC4 run (`tests/solve/adapt_indicators_test.cpp`, 4 cases, 342 assertions): the four-node fixture's edge values equal
the hand-computed ones (potential, both densities, doping by AR1, generation share 4 on node 3's edges); on the abrupt
1e17 diode the doping indicator alone marks only the nodes beside the junction and the potential indicator only nodes
within W/2 + h of it (W the abrupt-junction depletion width), zero excess beyond W/2 + 2h; on a MOS-C (p 1e17, 1.5 nm
oxide, n+ poly, 1 V, the quasi-static `equilibrium_poisson` solve: drift-diffusion stalls above threshold on a MOS-C,
`bias.hpp`) the electron-density indicator alone marks only nodes within one spacing of the inversion layer (n > N_A),
the surface node included. Fixture note: a constant profile's window along its depth axis (x by default) is ignored
(`doping_profiles.hpp`), so the x-windows take `axis=y`.

*Slice 3 — transfer, the loop, verification and the certificate (AD1, AD5, AD6, AD9, AD11, AD13; AC3, AC6, AC20
part).* `mesh/locate.hpp` (`PointLocator`: simplex meshes in 2D and 3D by orientation predicates through a bucket
grid, tensor meshes by `TensorCells::stencil`; outside: the nearest node); `GeneratedMesh::section` (a swept mesh's
cross-section simplices, node p of line j being j·n2 + p; not saved in files, so no file or digest changes);
`solve/adapt_transfer.hpp` (`transfer_state`: ψ, E_Fn, E_Fp and the density offsets c_n = ln n − (ψ + E_Fn)/V_T,
c_p = ln p + (ψ + E_Fp)/V_T interpolated over the old nodes of the new node's region; n, p rebuilt from them, so a
linear ψ and linear quasi-Fermi levels transfer exactly where c is constant, one material without narrowing);
`solve/adapt.hpp` (`adapt`, `next_refinement`, the builder callback `DeviceBuilder` and `MeshRefinement`). Without
adjoint weighting (slice 4) a goal's estimate is its last change between passes.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR2 | What the loop does when AD11 rejects a proposed stop. Found on AC6's first run: the reverse current's change between passes fell below tolerance at 70 nodes while its error against the 2× reference was 0.27%; AD11 rejected the stop (its verified error above the allowance; the estimate, 1.0e-12, within it), but the slice-3 loop as first written ended there, not certified. Continuing with indicator marks did not help (12 passes, the value settled 0.26% off: the indicators kept marking the same few nodes, so the missing accuracy lies where they see none; where exactly was not measured). | (a) End, not certified. (b) **Continue, the next pass marked by AD3 on the verified error's distribution** (AD12 c's wording, brought forward from slice 5 into the basic loop): per node of M, the largest over the adaptation points of \|ψ′ − ψ\|/V_T, \|ln n′ − ln n\|, \|ln p′ − ln p\| with M′'s solution interpolated at the node; uniform refinement (M′) when that marks nothing. (c) Continue with M′ itself. | **(b)**, ruled as recommended: it is AD12 (c) as frozen, reads solves only (OA6), and refines where the error is rather than everywhere; (c) multiplies 2D meshes by four per rejection. Slice 5 adds AD12's switch from the intelligent policy to this one; the basic loop has it from slice 3. No criterion changed. |

Runs (Release, `build/r3`):
- AC3: linear ψ and quasi-Fermi levels exact to 1e-12 tensor → tensor (5×4 to 9×7, not nested) and generated →
  generated (0.3 µm to 0.11 µm); the locator reproduces linear fields on tensor, 2D and swept 3D meshes. The 1e17
  diode at 0.5 V (41 → 97 nodes along x): currents equal to 1e-9, Newton 10 from scratch, 5 transferred; a lateral
  n+ p n+ MOSFET (gate 1 V, drain 0.5 V; 41×16 → 73×29): equal to 1e-9, Newton 12 from scratch, 8 transferred.
- AC6 (the legacy diode, weighted current on: the reverse current, 3.1e-9 A/cm², lies below the terminal current's
  resolution of about 4e-8): from 21 uniform nodes, 10 passes, AD11 rejecting two proposed stops (at 70 and 107
  nodes; AR2 then marked 30 and 64 nodes); stopped by the goal rule at **180 nodes** (the
  hand-made mesh: 250); I(−1 V) −3.09402e-9 against the 2× reference's −3.09688e-9 (9.2e-4, within AC6's 1e-3, a thin
  margin), I(0.6 V) 0.606567 against 0.606549 (2.9e-5); both certified; verification on 370 nodes in 0.02 s.
- AC20 part: a `max_nodes` stop and a `max_passes` stop are unverified (no verification run, verified error NaN);
  `verify = off` gives unverified after a goal stop. AD9: two runs bit-identical (node counts, goals, points,
  verified error, status). AC10 part: the refusals named in `adapt.hpp`.

*Slice 4 — the transposed solve, the goal adjoint, adjoint weighting and the estimate (AQ4, AD3, AD6; AC7).*
`linalg::LinearSolver::solve_transposed` and `can_transpose` (the equilibrated system transposed, (R A C)ᵀ y = C b,
x = R y, refined against the transposed componentwise backward error; SparseLU through Eigen's `transpose()` view,
PARDISO with iparm[11] = 2; gmres_ilut factorizes its sparse_lu fallback and reports `fell_back`, and without a
fallback cannot transpose, so `weighting = adjoint` is refused for it); `solve/adapt_adjoint.hpp` (`check_adjoint`,
`point_adjoint`, `goal_at`, `estimate_goals`); `TransferMode::exponential_fit` in `transfer_state`; the loop's
adjoint branch in `adapt.cpp`. The goal is the terminal current, or the weighted current G = T_c + Σ w_i B_i (Unit
32); its gradient is the small-signal currents' rows at s = 0 plus the weighted continuity rows, and the adjoint z
solves Jᵀ z = ∂G/∂u at the converged state.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR3 | How AD3's weights and AD6's estimate are computed from the adjoint (the design fixes the dual-weighted residual, not its discretization). Measured on the legacy diode (41 to 321 uniform nodes, against the solve on M′): M's adjoint interpolated onto M′ with M's state transferred by quasi-Fermi levels gave estimates off by orders of magnitude (the transfer's residual on M′ dominated); the exponential-fit transfer with M's adjoint still missed by factors of 5 to 50 at 0.6 V; a Newton-step (nonlinear) evaluation was worse. | (a) M's adjoint interpolated to M′. (b) **M′'s own adjoint at M's state prolongated to M′ by the Scharfetter-Gummel exponential fit** (n = Σλ_j n_j e^{(ψ−ψ_j)/2V_T} / Σλ_j e^{(ψ_j−ψ)/2V_T}, holes with the signs reversed; M's recombination densities N_s kept), one factorization of M′'s Jacobian per point: D = G′(u) − G_M − z′ᵀR′(u), the estimate 2\|D\| (p = 1, as AD11 without an observed order), the weights z′_q R′_q summed onto the coarse node with the largest interpolation weight in the same region and normalized over the semiconductor nodes, ω the largest over goals and points. (c) Uniform M′ solves each pass (AD11's cost every pass). | **(b)**, ruled as recommended: it is the dual-weighted residual on M′ with one linear solve instead of (c)'s nonlinear one, and the only form measured to work: at 0.6 V within 30% of G_M′ − G_M on every mesh; at −1 V within a factor 1.6 up to 161 nodes and 10 at 321 (a change of 1e-3 of the current); at 0.3 V the uniform meshes are not asymptotic (the change flips sign from 41 to 161 nodes), within 1.3 at 321. The cost is M′'s assembly and one factorization per adaptation point and pass (AC7: 2993 to 9570 nodes, 0.07 to 0.1 s). AD11 reuses the pass's M′. AR2's difference marking stays unweighted. No criterion changed. |
| AR4 | AC7's "fine uniform reference". The test first refined the start 4× (10449 nodes), then extrapolated 2×/4×/8× (5.18618, order 1.415), against which the adjoint estimate missed AC7's factor 3 (5.8). A 16× solve (164673 nodes, 12 minutes) moved the extrapolation to 5.20577 (4×/8×/16×, order 1.165): the uniform sequence is not asymptotic, so the reference is uncertain by about 0.4%, the size of the errors measured. | (a) Solve 16× in the test (12 minutes Release, far longer Debug). (b) **Record the 4×/8×/16× extrapolation as a measured constant, with the 2× solve repeated as the guard that the deck and physics are unchanged (to 1e-5).** (c) Keep 2×/4×/8×. | **(b)**, ruled as recommended: the finest reference measured, at suite cost; (c) uses a reference known to be 0.4% coarser. The reference's own uncertainty is reported to the owner: against (c) the estimate check fails, against (b) it passes; neither 32× nor any other independent reference was run. No criterion changed. |

Runs (Release, `build/r3`):
- linalg: transposed solves match factorizing Aᵀ (real and complex, transpose not conjugate; backward error within
  1e-9 measured 2.4e-10); PARDISO transposed with MKL; gmres_ilut falls back, or refuses without a fallback.
- AQ4: the goal's gradient against central differences on the diode (terminal and weighted current at −1, 0.3 and
  0.6 V): within 1e-6 of the terms' magnitude in five cases; **the weighted current at −1 V differs by 5.8e-3** of the
  derivative (open finding: not isolated; not SRH, Auger, radiative or narrowing alone; the test bounds it at 1e-2
  in that case only).
- AC7 (a MOSFET in a p-well, 1 µm × 1.5 µm, gate 1.5 V, drain 0.5 V; I_D's reference 5.20577, AR4): from 693 nodes,
  feature-only stopped by the goal rule at 6834 nodes (I_D 5.18975, 0.31%, certified); adjoint-weighted at **2214
  nodes** (I_D 5.18151, 0.47%, certified; estimate 0.0271, AD11's verified error 0.0156, actual 0.0243, estimate /
  actual 1.12); node ratio 0.32 (AC7: at most 0.75); y-lines within 0.15 µm of the well junction 3 (the initial
  3) against feature-only's 19. 29 s.
- Open finding for the owner: AD11's verified error \|G_M′ − G_M\| / (2^p − 1) estimates M′'s error, not M's (M's
  is 2^p times it, \|G_M′ − G_M\| 2^p / (2^p − 1)); frozen, not changed.

*Slice 5a — the decision layer's prediction, guard and fallback (ID1–ID3, ID7 part, AD12; AC14, AC17, AC21).*
`solve/adapt_policy.hpp` (`Policy`, `predict_sizes`, `edge_exceeds`, `choose_points`, `PolicyOptions`); the loop's
`policy = intelligent` branch in `adapt.cpp` (the coarse sweep and point choice, the prediction's rounds before pass 0,
the miss rule, AD12's three triggers, `MeshSource` per pass, the fallback and misses in the result). The default stays
`policy = basic` (ID9: not the default until AC14–AC17 pass). Slice 5b adds ID4–ID6 and ID5's cost model.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR5 | The constructions the decision layer's rows leave open. | (a) as below; (b) other thresholds. | **(a)**, ruled as recommended; no criterion changed. **Split:** slice 5 is two commits, 5a (ID1–ID3, AD12, ID7's per-pass record) and 5b (ID4–ID6, ID5, ID7's candidates). **Miss (ID1):** a proposal "lowers the estimate" when it falls by at least the miss margin (10%) of the previous pass's; the estimate compared is the largest over goals and points of estimate / allowance (the adjoint's, or the last change without weighting); a proposal that refines nothing is a miss without a solve; pass 0's mesh is judged only against ID3's coarse sweep's adjoint estimate. **Prediction (ID2):** applied before pass 0 in rounds (each round halves sizes toward the targets, at most 12, recomputed on each round's device), stopping before ID6's growth cap (2.5× the start's DOFs; the prediction is pass 0's proposal) or max_nodes; a round that would cross the cap is applied in part, to the edges of largest length / target ratio, the share found by bisection on the sorted ratios (at most 12 trial meshes), and the rounds then stop. Targets are anisotropic on tensor and swept axes (an edge's length projected on the feature's normal) and isotropic in 2D size fields. Under gates the depletion region is added (L_D / 2 within W_d, the junction's own rule; the inversion length alone covered almost nothing the basic loop refines on the MOS-C). Impact ionization and band-to-band tunnelling add nothing beyond their junction's target. **Points (ID3):** a point where a goal is at or below its current resolution is not eligible (its allowance tol·r cannot be met, AD6); the coarse sweep is pass 0's solve when the prediction refined nothing. **Hand-written prediction:** `PolicyOptions::prediction` replaces ID2 and every later proposal (AC17's wrong field). **AC17 on the diode** checks the guard; its reference is the 4× hand mesh (the 2× moves 4.7e-4 from it at −0.5 V, half the tolerance); its time bound is measured on the 2D MOSFET, the device whose substrate the criterion names. **AC13:** AC7 (28 minutes in Debug) and AC14 are hidden gates run in the Release configuration only, as the solve layer's slow gates are. |
| AR6 | AC14's devices and its measured miss. | (a) as below; (b) the generated 2D MOSFET. | **(a)**, ruled as recommended; no criterion changed. **Devices:** the diode (AC6's run), the MOS-C (substrate current, degenerate zero goal), and the 2D MOSFET on a 25 × 11 tensor mesh (I_D at Vd 0.05 V, Vg 0.5 and 1 V). The generated (isotropic) 2D MOSFET was measured with both policies to 12 passes: about 85 k nodes, I_D still moving, 47 minutes, not converged (as Unit 14's runs); it stays a probe (`[.adapt_probe]`) and is a risk for AC5 (slice 7). **Metric:** recall and precision count node–axis pairs (a node refined along x and y is two); on generated meshes one isotropic axis, refined below 0.9 of the start size. **Measured 2026-11-12** (slice 5a, before ID4–ID6): recall 0.87 / 0.83 / 0.51 (≥ 0.8 fails on the MOSFET: the prediction refines across its features, the basic loop both axes; per node it was 0.94); precision 1 / 1 / 0.996; false-positive DOFs 0 / 0 / 102 of 26460; DOFs 1116 vs 540 and 504 vs 216 (≤ 1.1× basic fails on the diode and the MOS-C), 26460 vs 31752; passes 11 vs 10, 2 vs 2, 5 vs 6 (saved 0, fails); MOSFET 76 s vs 100 s. AC14 is marked `[!shouldfail]`, re-measured after slice 5b; owner ruling then. |

*Slice 5b — benefit, cost model, online fits and allocation (ID4–ID6, ID7's candidates; AC15, AC16, AC18).* The
factor fill nnz(L + U) is reported by every backend (`FactorizationReport::fill`: SparseLU's nnzL + nnzU, PARDISO's
iparm[17], ILUT's factors) and kept per Newton solve with the Jacobian's nonzeros (`ConvergenceRecord`, not saved or
observed). `solve/adapt_policy.hpp` adds `fit_cost_model`, `fit_orders`, `order_of` and the budgets; `adapt.cpp`
the candidates, their benefits, work and predicted seconds, the allocation, the realized ratios and the refits;
`AdaptPass` the per-pass DOFs, nonzeros, fill, fitted models, predicted and listed candidates.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR7 | The constructions ID4–ID6 leave open, and a conflict between ID4 and AD9. | (a) as below; (b) rank by the timing-fitted seconds (not bit-identical, AD9 and AC18 fail by construction). | **(a)**, ruled as recommended; no criterion changed. **Determinism (owner to confirm):** ID4's cost is seconds from a timing-fitted model, so ranking by it would make the mesh depend on the machine's timings, against AD9 and AC10/AC18. The allocation ranks by a deterministic work in ID5's form at its default exponent, N_newton ((W + dW)^1.5 − W^1.5), W = nnz + fill (counted, not timed), dW the candidate's added nonzeros at the mesh's mean degree and the candidate's unknowns per node, with their fill at the current ratio; the timing-fitted model predicts each candidate's and each pass's seconds (reported, gated by AC15) and enforces only a time budget (wall-clock by nature). **Adjoint only:** ID4's benefit is "in the goal's units through the adjoint"; without weighting the estimate is the last change, which a large refinement raises (measured: the tensor MOSFET's proposals took it from 2.5 to 3.5 at 19 k nodes, then the fallback refined to 51 k nodes against basic's 10.6 k), so the allocation runs with `weighting = adjoint` only, and without it the layer is slice 5a's (ID2, ID3, AD3's marks). **Candidates:** an axis interval (tensor axes, the swept axis; one line added per halving) or a node (a generated 2D size field or a swept cross-section, its size halved; 3 nodes added in 2D, 3 per line swept). An interval can be split to three levels in one pass, each level a greedy entry of benefit E (S_c / S) 2^−(l−1)p (1 − 2^−p) and 2^(l−1) times the first level's nodes (without levels every candidate with an excess was taken while the target is unreachable, θ → 0: the tensor MOSFET then needed 2× basic's DOFs); a node one level (a size field is halved at most per build). S_c is the candidate's weighted excess (ω η: intervals over their edges along the axis, with ω the larger end's; nodes their own; a swept cross-section node over its lines), S the sum over candidates, E the relative estimate (AR5). **Target and caps:** take entries in decreasing benefit over work (ties by construction order; a level after the one before it) until the predicted estimate ≤ min(1/2, E/2) (tol/2 at most; at least halving, so a proposal is not a miss by construction), the DOF growth cap (2.5× per pass) or the budget. Both caps hold on the built mesh: grading adds lines and sizes the prediction misses (measured 1.5× on the well MOSFET), so a proposal built over a cap is allocated again with the measured inflation (built over predicted added DOFs), at most three times, and the inflation is kept for the next pass. Budget mode (`dof_budget`, `seconds_budget`): the same order until the budget is spent; a pass that fits nothing, or a built mesh over the budget, stops the run (`max_nodes`, unverified); the prediction's rounds (ID2) also stop at the budget. **Realized ratios (ID5):** a taken candidate's local level, the ω-weighted mean over its place (an interval's slab; a node's disc of its local size, in its region; a swept cross-section node's cylinder) of the node level max over edges of raw / threshold for its indicator kind (the excess vanishes below the threshold, the level does not), on the next pass over this one; per halving its l-th root. Measured for the first 200 taken per pass, the first 20 listed (ID7). **Fits:** the order per (kind, region), the mean of −log2 ratio with a ridge of one sample toward 1, clamped to [0.5, 3]; the cost model per Newton iteration a nnz^b + c fill + d DOF with mesh and transfer per node, one sample a alone at b = 1.5, two a and b, three or more b on a grid of 0.05 in [1, 2] with a, c, d non-negative least squares of the relative error and a ridge (0.01 / n)(b − 1.5)² (0.05 fixed outweighed the misfit on exact data); a pass's time predicted from the earlier passes with the last pass's Newton count at this mesh's measured nonzeros and fill. The node-only model T = αN^β (AC15's comparison) is fitted beside it. **AC16 devices:** AC5's run is the generated MOSFET's (AR6), so the tensor MOSFET (AC14's) stands in, adjoint-weighted; "the same verified accuracy" is both certified at the tolerance, the errors against the 8× uniform reference reported (4× and 8× differ by 2.8% at Vg 0.5 V: the reference is not converged to 1%, and both certified results are 3% from it; for AC20, slice 7). Basic at a budget is basic with `max_nodes` at the budget's nodes. **AC15 devices:** the tensor MOSFET for AC5's run, the adjoint diode for AC6's; AC8 (the FinFET) is slice 7's. **AC18:** "bit-identical, report included" excludes the measured seconds and what is fitted to them; the report's numbers are recomputed from the saved passes (orders, benefits and realized reductions, the predicted estimate, both fitted models and the predicted seconds). |

**Slice 5b measured (2026-11-13, Release, r3):**
- AC16 tensor MOSFET (adjoint, tol 1e-2, I_D at Vg 0.5 and 1 V): basic 8 passes, 15 732 DOFs, 23 s; intelligent 4
  passes, 9 270 DOFs (0.59×, ≤ 0.8 holds), 9 s; both certified, errors against 8× 3.1% and 3.4%. Budget 7 866 DOFs:
  intelligent 7 614 DOFs, error 2.8945%; basic 3 960 DOFs, 2.8969% (holds, barely). **Passes.**
- AC16 well MOSFET (AC7's deck, adjoint, I_ref = 5.20577): basic 6 passes, 6 642 DOFs, error 0.47%; intelligent 3
  passes, 12 744 DOFs (1.92×), error 0.32%, relative estimate 0.013 at the stop against basic's 0.52: the first
  allocation, growth-capped from E = 9.5, overshoots the accuracy. Budget 3 321 DOFs: intelligent 3 255 DOFs, 6.5%;
  basic 3 240 DOFs, 7.2% (holds). **Fails in tolerance mode**; `[!shouldfail]`, owner ruling.
- AC15: rank correlation 0.97 (MOSFET, 14 candidates) and 0.997 (diode, 37) (≥ 0.6 holds); the cost model from
  the third pass within 30% (worst 5% MOSFET, 15% diode; holds) but not better than node-only (mean 3.3% vs 2.7%,
  8.0% vs 7.6%: both small and within timing noise on 2 and 12 passes; fails); the adjoint diode's orders 1.10
  (potential, 32 samples), 1.09 (electrons, 2), 0.81 (doping, 4) and 1.64 (holes, 2; outside [0.8, 1.3], fails).
  `[!shouldfail]`, owner ruling. The adjoint diode does not converge in 14 passes with either policy (its estimate
  moves between 2 and 15 on meshes a few nodes apart; basic adjoint likewise: 571, 85, 105, 2.6, 25, …); that is the
  weighted-current gradient finding of slice 4, not the allocation.
- AC18: two runs bit-identical (meshes, goals, estimates, candidates, fits, certificates), the report recomputed
  from the saved passes. **Passes.**
- AC14 is unchanged by 5b (its runs are feature-only); with adjoint weighting the intelligent loop takes 4 passes
  against basic's 8 on the tensor MOSFET (AC16's run).

**Owner rulings after slice 5b (2026-11-13; owner, "Decisions for you: choose recommended").** No criterion is
loosened or re-worded. (1) AR7's determinism construction is confirmed: the allocation ranks by the counted work
N_newton ((W + dW)^1.5 − W^1.5), and the timing-fitted model only reports (AC15) and enforces the time budget.
(2) The allocation runs with adjoint weighting only; feature-only runs keep slice 5a's layer. (3) AC14 (AR6), AC15
and the well-MOSFET AC16 stay `[!shouldfail]` with their measurements above. They are re-measured once slice 7's gates
are in, and any still unmet at the unit's close is recorded as a limitation in ARCHITECTURE.md 11.3, as for Units
46, 48, 53 and 54. (4) AR4: AC7 keeps the recorded 4×/8×/16× reference 5.20577. Its uncertainty (about 0.4%) is stated
beside AC7 and re-examined with the tensor MOSFET's unconverged 8× reference under AC20 in slice 7.

**Slice 6 (AD7, AD8, AD13; AC9, AC10).** `[adapt]` in the deck (deck/build.hpp), the runner's adaptive path
(deck/run.cpp), `results::AdaptReport` (results/adapt.hpp) filled by `solve::adapt_report`, its file
(io/adapt_report.hpp, content kind 8, `.nitadapt`), the identity stamp (`solve::stamp_adapt`), the certificates in
`nitcad run`'s output and the CSV table, exit status 3, and the app's deck vocabulary. 6b adds `nitcad plot`'s and the
app's views of the passes.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR8 | The constructions AD7, AD8 and AD13 leave open. | (a) as below; (b) implement the derived goals now. | **(a)**, ruled as recommended; no criterion changed. **Goals:** `goal = current <contact>` only; `vt`, `gm` and `swing` are refused with their line ("not adapted in this unit"). The loop, its adjoint (AR3) and AD11 are built on terminal currents, and a derived goal needs its own adjoint (a functional of several points) and resolution rule: a later unit, **confirmed by the owner (2026-11-13)**. `predictor` (17b) is refused. **Defaults:** `weighting = none`, `policy = basic` (ID9: intelligent is not the default until AC14–AC17 pass, and AC14 and AC15 do not), `verify = richardson`, `certify = report`, `tolerance = 1e-3`, `max_passes = 6`; `at = auto` (absent: the same): the sweep's first and last points, or, with the intelligent policy, up to `max_points` (default 3) chosen by ID3; `at = <indices>` counts the sweep's points from 0. `mode = budget` needs the intelligent policy and `budget_dofs` and/or `budget_seconds` (the basic loop takes `max_nodes`); the budget keys are refused in tolerance mode. **Scope refusals (AC10), each with its line:** an analysis other than `[sweep]` (AQ3: one mesh per sweep), a restart (the checkpoint is the deck's mesh's), a generate3d, file or unswept extrude3d mesh (AQ2), `weighting = adjoint` with `check_adjoint`'s refusals (a backend that cannot transpose and no fallback, electrothermal, energy balance, nonlocal tunnelling) or a non-ohmic goal, a goal on a gate or electrode, a second goal on one contact, the adaptation's outputs without `[adapt]`. **Rebuild (AQ1):** the builder is `deck::build(deck, refinement)`, the whole deck rebuilt per pass: a tensor axis the loop refined is taken whole (its straddles already in it; their spans protected), a generate2d mesh and a swept cross-section take the size field, a swept mesh the refined lines; an empty refinement gives the deck's mesh (AD10). **Outputs (AD8):** `sweep`, `device`, `mesh`, `checkpoint` and `table` are the final mesh's; `adapt = <file>` the report; `adapt_meshes = <directory>` each pass's device and mesh (`pass_<k>.nitdev`, `pass_<k>.nitmesh`, written as each pass is solved through `AdaptOptions::on_mesh`): the passes' meshes the app and `nitcad plot` show. The report holds the run record, the settings, per pass the mesh digest, source, goals, estimates, marks by indicator and region, Newton counts, nonzeros, fill, the decision layer's predictions, candidates, orders and cost model, the times; the certificates; the stop, fallback and verification; its status is the worst certificate's. **Identity (AD8):** the final sweep's `input_identity` continued by FNV-1a 64 over "adapt", each setting's name and value (bit pattern) and the final mesh's digest, the settings appended as `adapt.*` (`solve::adapt_settings`); a run without `[adapt]` is unchanged (AD10). **Certification (AD13):** `nitcad run` prints each pass and each certificate, with "an estimate and an independent Richardson check, not a proof"; the CSV table gains per goal `<c> certificate`, `<c> estimate A`, `<c> verified error A`, `<c> allowed A` and `<c> order` after `substeps` (filled at the adaptation points); with `certify = required` a result that is not certified exits **3** (2 stays an early stop), its outputs written. |

**Slice 6 finding (not changed):** the deck tests' `diode_device` fixture (tests/deck/fixtures.hpp) is not a diode. Its
profiles `acceptor constant 1e17 x=-1:1um` and `donor constant 1e17 x=1um:3um` window the depth axis x, which a profile
ignores (doping_profiles.hpp: "the depth axis's entry is ignored"), so both dopants cover the whole device and it is a
compensated resistor (the anode current is linear in the bias, 1.46 mA at 0.2 V). Its tests compare a deck against the
same hand-built device, so they still hold. Slice 6's deck uses a lateral window (`axis=y`) on a 21 × 2 mesh, which
gives a junction (2.6 nA at 0.4 V). Owner to decide whether to correct the fixture (every deck test that uses it would
change its numbers).

**Slice 6b (AD8: the views).** `nitcad plot --view adapt` draws each pass's estimate over its allowance against its
DOFs (log-log, the tolerance the line at 1), from the deck's `[output] adapt` or a `.nitadapt` file; `--view passes`
tiles each pass's mesh from `adapt_meshes` (`--columns`, 4; `--size` per tile, 640×480), each pass file checked against
the report's mesh digest. Every other view of a deck with `[adapt]` draws the **final** mesh, read from its `device`
(and `mesh`) output: the results were solved on it, not on the deck's start mesh, so `plot` (and the app) would
otherwise refuse them by digest. The pass files are named by `deck::pass_device_file` / `pass_mesh_file`. The app:
`Workspace` carries the report, the passes' directory and, for a pass shown alone, its index; `workspace_from_run`
takes the final mesh's device from the run, `workspace_from_plan` from the device output or the last pass's file
(whichever has the report's final digest; an error when neither does: the report's results cannot be paired);
`workspace_of_pass` loads one pass's mesh without results; the run record gains the adaptation (status and stop, "an
estimate and an independent check, not a proof", passes and nodes, policy and weighting, fallback, verification, each
certificate) and the pass shown. The results page lists the passes (pass, mesh source, nodes, DOFs, estimate over
allowance, marked, seconds, notes) for an adaptive run only; activating a row shows that pass's mesh in the device
view with the results kept aside, and the last pass brings them back (`MainWindow::results_show_pass`). Tests: the
plot views, their refusals and a changed pass file (tests/deck/adapt_test.cpp), the workspace (tests/app/
workspace_test.cpp) and the window (tests/app/analysis_window_test.cpp). Slice 6 is committed as one commit (6a and
6b), its suites run on both together.

**Slice 7 (gates and cost: AC5, AC8, AC11, AC20, AC12; ARCHITECTURE 6.12.4 and 11.3).**

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR9 | The tensor MOSFET's reference (AC16's comparison, AC20's on AC5's stand-in, AR6). The 8× uniform reference of slice 5b is not converged. | (a) per-axis ladders, each axis extrapolated, recorded with a 2× guard; (b) a finer uniform ladder (32×: 590 k nodes, hours); (c) keep 8×. | **(a)**, ruled as recommended; no criterion changed. **Measured (Release, 2026-11-13), I_D (A/cm) at Vg 0.5 and 1 V.** Uniform 1×–16× (275–61 985 nodes): 0.13051, 0.13257, 0.13214, 0.12846, 0.12575 and 0.28192, 0.30445, 0.31875, 0.31947, 0.31763: not monotone, because the y error (the inversion layer: 16× is 1.25 nm) and the x error have opposite signs and different orders. Per axis: at kx = 4, ky = 32, 64, 128: 0.122849, 0.122425, 0.122309 and 0.305497, 0.304497, 0.304174 (orders 1.87, 1.63); at ky = 32, kx = 4, 8, 16, 32 (to 246 849 nodes, 14 min): 0.122849, 0.124059, 0.124385, 0.124483 and 0.305497, 0.312421, 0.314917, 0.316039 (orders 1.89 then 1.73; 1.47 then 1.15). The y step 32 → 64 at kx = 8 is within 4% (0.5 V) and 12% (1 V) of that at kx = 4: the errors add. Each axis extrapolated at its last observed order, the y correction scaled to the finer x: **0.12393 (uncertain by about 0.1%) and 0.3153 (about 0.3%: the x order still drifts)**. The 8× values are 3.7% and 1.3% high, which is slice 5b's "certified runs 3% off". The values are recorded in adapt_intelligent_test.cpp (`mosfet_ref`), with the 2× solve (0.1325677205, 0.3044490795) as the guard that the device and the physics are unchanged (as AR4). AC16 now reports its errors against them instead of solving 8× and 4× each run. |

**AC11 (measured 2026-11-13, Release, another 246 k-node solve running beside it; timings re-measured idle at
close).** The estimate's factorizations of J′ on M′ and its largest transposed solve are timed (`PointAdjoint`,
`GoalEstimate`, `AdaptPass::seconds_estimate_factorize` and `seconds_transposed`, the report row and `.nitadapt`; both
excluded from the AD9 comparison as measured seconds). The gate compares one transposed solve against one
factorization, a lower bound on a Newton iteration's linear solve (which also solves). Tensor MOSFET, adjoint: basic 8
passes (275 → 5 244 nodes; nnz 7 615 → 150 476; fill 49 890 → 2 671 164; 18–22 Newton iterations per pass), the
transposed solve 2.2–5.4% of a factorization (0.047 s against 1.32 s on the last pass); verification 21 436 nodes,
14.9 s against the final pass' 1.6 s. Intelligent 4 passes (663 → 3 090 nodes), the transposed solve 2.5–5.2%;
verification 12 900 nodes, 5.6 s. Adjoint diode, 14 passes (21 → 143 nodes), 3–13%; not verified (the max_passes stop
of slice 5b). **Passes** (the AC5 and AC8 rows with their own runs below).

**AC20 on AC5's stand-in (measured 2026-11-13; AR9's reference).** Tolerance 1e-2. Basic: certified at both points,
actual errors 0.156% and 0.431%; the verified error 0.33 and 0.87 times the actual (p 1.86 and 0.5), the adjoint
estimate 1.65 and 1.17 times. Intelligent: certified, 0.452% and 0.198%; verified 1.88 and **0.038** times the actual
(p 0.5 and 2.18), estimate 1.59 and 0.38 times. Each run stopped by a `max_nodes` budget one node above its first pass
is unverified at both points (errors 5.3–11%). **Passes.** Finding (not a failure, recorded): AD11's verified error
can underestimate the actual error by a large factor when its observed order is high (2.18 here, where the order
clamp [0.5, 3] lets it extrapolate the last change down by 2^2.18 − 1 = 3.5 and the change itself is small); the
certificate held only because the actual error was a fifth of the allowance. AC7's case now also checks that a
certified result is within its tolerance of AR4's reference.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR10 | AC5's reference and its "uniform-refinement ladder needs" (AC20's "AC5's 4× uniform"). | (a) AR9's converged reference (the same continuum device), the 4× generated value reported beside it; the ladder's nodes bounded below from its h0/4 and h0/8 runs at an assumed second order; (b) as frozen: the generated mesh's own 4× reference, the ladder measured to 1%. | **(a)**, ruled as recommended; **confirmed by the owner (2026-11-13)** (it changes the reference AC5 and AC20 name, toward a more accurate one, not a looser one). The generated MOSFET is the tensor MOSFET's device (geometry, doping, contacts and gate identical; only the mesh differs), so its mesh-converged I_D is AR9's. Its 4× uniform value (5 nm, 20 295 nodes) is −13.8% and −19.3% from it (8×, 80 660 nodes: −6.1%, −11.0%), so as a reference it would certify nothing meaningful at 1%, and the ladder cannot be measured to 1% (it converges below first order between 4× and 8×; second order from 8× on would still need about 890 k nodes; the observed order about 1 M × 30). Vg 0 (above its resolution, 4.8e-9 A/cm) gets the same per-axis construction: x at ky = 32 (kx 4, 8, 16): 7.08387, 7.08502, 7.08562e-5 (order 0.94); y at kx = 4 (ky 32, 64, 128): 7.08387, 7.03040, 7.01701e-5 (order 2.0), the y step at kx = 8 equal to 0.1%: **7.0150e-5 A/cm** (about 0.05%). |

**AC5 (measured 2026-11-13, Release; generated MOSFET from h0 = 20 nm, 1 323 nodes; I_D at Vg 0, 0.5 and 1 V,
tolerance 1e-2, adjoint-weighted; against AR9 and AR10).** Uniform: h0/2 5 147 nodes (−23.6%, −27.4% at 0.5 and 1 V),
h0/4 20 295 (+12.9% at 0 V, −13.8%, −19.3%), h0/8 80 660 (+2.7%, −6.1%, −11.0%). **Basic** (12 passes, 17 s): 1 323 →
2 888 nodes, errors −15.3% and −19.9% at the last pass (θ = 0.3 marks a few nodes; growth about 10% a pass), relative
estimate 22.9, max_passes, unverified. **Intelligent** (9 passes, 40 min): the prediction's 3 303 nodes, then
proposals 7 202, 14 823, 30 104, 68 898 (+0.35%, −0.36%, −1.33%), **98 593 (+0.37%, −0.20%, −0.89%: inside 1% at every
point)**, 120 824, 129 977, 146 084 (+0.30%, −0.19%, −0.83%); relative estimate 1.05 at the stop (the estimate 1.25 times
the actual error at 1 V), so it stopped at `max_nodes` (150 000), unverified at every point. **98 593 nodes against at
least 890 k for the ladder (0.11 of it; ≤ 1/3 holds) with the intelligent policy; the basic loop does not reach 1%
within 12 passes.** The gate (`[.ac5]`, manual like `[.uc8]`: about 40 minutes) repeats the h0/4 and h0/8 runs for the
bound and the intelligent run to 100 000 nodes, and checks AC20 on it.

**AC8 (measured 2026-11-13, Release, PARDISO on 6 threads; `[.ac8]` in tests/deck/complex_finfet_test.cpp).** The
start: decks/finfet_3d_swept.deck with its axis lines about 10 nm apart (0, 20, 25, 35–50 every 5, 60, 65, 85 nm, the
gate edges straddled as the deck) and its cross-section unstructured at up to 12 nm (the structured frame unchanged):
4 718 nodes. `[adapt]`: goal the drain current, tolerance 1e-2, adjoint-weighted, `at = auto` (0 and 0.8 V for basic,
ID3's points for intelligent), max_nodes 9 367 (UC8's budget). The coarse tensor mesh (9 367) and the 2× reference
(45 375) re-solved as UC8. **Basic:** 2 passes (4 718 → 7 678 nodes; the next would exceed the budget), relative
estimate 175 then 139, stop max_nodes, unverified; UC8's rule at **9 of 17** points (fails at 0.1, 0.15 and 0.35–0.6
V: on-current 3–5% low at 0.4–0.55 V, the tensor within 0.6%), V_T and the swing within their bounds, max g_m −8.3%
outside. **Intelligent:** the prediction's 7 350 nodes and no further pass, relative 125, unverified; **6 of 17**,
V_T and max g_m outside, the swing within. UC8's hand-tuned deck (Unit 54, 8 300 nodes) met 12 of 17. **The measured
target is missed** (it was stated as a target, not a gate). The gated parts hold: every point converges, the budget
holds, both runs unverified (AC20). The cause measured: on the swept 3D mesh the adjoint estimate is 84 times the
actual error at 0 V (2.98e-10 against 3.56e-12 A) and 9 times at 0.8 V (1.29e-6 against 1.37e-7 A), so the loop sees
itself 100× from its tolerance and the budget allows one or two passes, spent where the overestimated weights point.
AC11 on AC8: per pass 6 734 and 11 550 DOFs, 181 820 and 302 032 nonzeros, fill 6.46 M and 13.1 M, 109–115 Newton
iterations over the 17 points; solve 5.8 and 10.1 s, the estimate 3.1 and 5.2 s (factorizations 1.0 and 2.1 s,
transposed solve 0.024 and 0.041 s: 2% of a factorization); no verification run (budget stops).

**Re-measured after slice 7's gates (2026-11-13, Release, idle machine; the owner ruling after slice 5b).**
`[adapt_intelligent],[adapt_goal]`: AC7 passes (adjoint 2 214 nodes against feature-only 6 834, 0.32; error 0.47%
against AR4). AC11 and AC20 pass. **AC14** unchanged (its runs are feature-only): recall 0.87 / 0.83 / 0.51, precision
1 / 1 / 0.996, DOFs 1116 vs 540 and 504 vs 216, passes saved 0: **unmet**. **AC15**: rank correlation 0.97 and 0.997
(holds); the diode's orders 1.10, 1.09, 0.81 and **1.64** (holes, 2 samples), and the cost model not better than
node-only within timing noise: **unmet**. **AC16** tolerance mode on the tensor MOSFET: 9 270 against 15 732 DOFs
(0.59), both certified, errors against AR9 0.45% and 0.43%: **holds**. **AC16 budget mode on the tensor MOSFET, against
AR9's reference: intelligent 0.6531% (7 614 DOFs) against basic's 0.6507% (3 960 DOFs) at a 7 866-DOF budget —
unmet** (a tie within AR9's 0.1% uncertainty; slice 5b's pass was against the unconverged 8× reference, 2.8945%
against 2.8969%). The case is split: tolerance mode stays a gate; budget mode is `[!shouldfail]` with this
measurement, the criterion unchanged. **AC16 on the well MOSFET** tolerance mode 12 744 against 6 642 DOFs (1.92×):
**unmet**; its budget mode holds (6.5% against 7.2%). Per the ruling, the unmet criteria at the unit's close are
limitations in ARCHITECTURE.md 11.3: AC14, AC15, AC16 (well tolerance mode; tensor budget mode), and AC8's target.

**Slice 7 suites (2026-11-13, at ca22846):** Release 217/217 (980 s, `-j 6`); Debug 198/199, the one failure
`app.live_check` (its wake count, line 53: the known flake under `-j 6`, untouched by this unit), which passes alone
three times in a row (`--repeat until-fail:3`).

**AC12 (hand mutants, 2026-11-13).** Each mutant one exact-string change, the solve, deck and mesh libraries and
`nitcad_adapt_test`, `nitcad_deck_test` and `nitcad_mesh_test` rebuilt (Release, r3), then the fast suites: every
non-hidden adapt case, the deck's adapt cases and `[size_field]`; killed if a build or a test fails. **Round 1, 37
mutants:** 23 killed; 13 survived; 1 did not build (ID5's, an unused variable under `/WX`; replaced). The survivors were
rules the end-to-end runs never pinned: AD1 (solve from scratch, not the transferred state), AD4 (halving 0.5 → 0.7;
the `min_size` skip dropped), AD9 (the tie-break reversed), AD11 (the default order 1 → 3), AD12 (the miss count not
reset), AD13 (`&&` → `||` in the status; the allowance doubled), ID1 (the margin's sign), ID3 (one point, not
`max_points`), ID6 (the target E, not E/2; the growth cap ×10), ID7 (2 listed, not 20). **Closed by tests:** the rules
are functions in solve/adapt.hpp, applied by `adapt()` in place of the inline code (no behaviour change; the fast
suites' numbers unchanged): `allowed_error`, `richardson_error`, `certificate_status`, `refined_size`, `heavier_node`,
`proposal_missed`, `MissGuard`, `allocation_target`, each tested exactly in adapt_loop_test.cpp ("the loop's
rules"), with `CostModel::per_iteration`'s form; AC6 checks that every pass after the first was transferred onto
(AD1); AC17 that `max_points` points are chosen (ID3) and that the prediction's rounds and every proposal hold the
growth cap (ID6); AC18 that the first 20 taken are listed (ID7) and that an allocation away from the cap stops at
its target (ID6; on that diode E never falls below 1, where the target's two forms differ, hence the exact test).
**Round 2:** the 13 re-aimed at the rules' new places, and ID5's (the DOF term dropped): 12 killed; ID5's and ID6's
target survived; **round 3** after their tests: both killed. **Every design row has at least one killed mutant and
none survives:** AD1 (2), AD2 (2), AD3 (2), AD4 (3), AD5, AD6 (3), AD7, AD8 (2), AD9, AD10, AD11 (4: AC20's skipped
refinement, wrong p twice, wrong mesh compared), AD12 (3), AD13 (3: AC20's status logic twice, the allowance), ID1,
ID2, ID3, ID4, ID5 (2), ID6 (2), ID7, ID9. **AC12 met; AC20's mutant clause met.** The runner first restored each file
with `git checkout`, which discarded uncommitted edits to the mutated file in round 2's first attempt; it now restores
from a copy, and each round runs on a committed tree.


---

### Unit 17b — cross-run learning for adaptive refinement (started 2026-11-13; on `mesh/cross-run-learning`)

**Start (owner, 2026-11-13: "start Unit 17b, cross-run learning").** The design is Unit 17's frozen rows ID8, ID9
(17b's part), ID10, AC12 (ID8, ID10), AC19 and AC21's predictor case, with amendment OA5; nothing there is changed
here. What follows are the constructions those rows leave open, each put as a question with a recommendation and
ruled as recommended under the owner's standing instruction (2026-11-11, "always choose recommended"). None loosens a
criterion; LB10 narrows what one design row stores and was confirmed by the owner.

| # | Open construction | Options | Ruling |
|---|---|---|---|
| LB1 | The base. The split says "a separate branch and unit after 17 merges"; PR #81 (Unit 17) is open. | (a) `mesh/cross-run-learning` from `mesh/adaptive-refinement` at 6ffad0b (stacked); its PR to main once #81 merges (or retargeted then); (b) wait for #81. | **(a)**, ruled as recommended. Unit 17 is accepted and pushed; 17b reads its code, not main's. Any later change to #81 is merged into this branch. |
| LB2 | The sample and its label ("node features → realized benefit", ID8). The predictor's output is a first mesh, i.e. sizes. | (a) **One sample per edge** of every solved pass's mesh; its label the realized refinement ℓ = max(0, log2(L_e / h_f)), L_e the edge's length and h_f the final mesh's local size along the edge's direction at its midpoint (the extent, projected on the edge's unit direction, of the final cell containing the midpoint, by Unit 17's `PointLocator`); an edge whose midpoint lies outside every final cell is skipped. The realized benefit is the refinement the converged adaptation found worth its nodes there. (b) Per refined candidate, the adjoint estimate's realized reduction (ID5's samples): defined only where something was refined and only with adjoint weighting, and not a size. | **(a)**, ruled as recommended. ℓ is defined on every edge, with any weighting, and maps straight to a size (LB4); the per-candidate reductions stay ID5's online fit. "Benefit prediction" in AC19 and ID10 is the prediction of ℓ. |
| LB3 | The features (ID8's list) and the schema. | (a) **Schema 1**, per edge from its end nodes (the mean of the two, or the edge's own): region kind (semiconductor, insulator: a category); material (its name's index in the model's sorted material names: a category); log10 max(\|N_D − N_A\|, n_i); \|∇ log10 max(\|N\|, n_i)\| along the edge (per µm); log10 of the distance to the nearest junction point, insulator-interface point and contact node (cm; 1 cm when there is none); log10 of ID2's W and L_D of the nearest junction point and of λ of the nearest interface point (`predict_sizes`' construction; 1 cm when none); the analysis' largest junction and gate bias differences (V, ID2's); the goal kind (0: a terminal current, AR8); log10 L_e; \|cos\| between the edge and the nearest feature's normal; the nearest junction's distance over its W and the nearest interface's over its λ. 17 features. (b) ID8's list alone. | **(a)**, ruled as recommended. The additions are what the label depends on that ID8 did not list: the edge's own length and direction (an edge across a junction and one along it are refined differently, Unit 17's projected rule), and the two scaled distances (a device family varies L_D and λ; trees cannot form ratios). A change to any feature's definition is a new schema version (ID10). |
| LB4 | From a prediction to the first mesh. | (a) Per edge the target L_e 2^−ℓ̂ where ℓ̂ ≥ 0.5 (none below); per node the smallest over its edges, with that edge's direction as the normal (Unit 17's `SizeTargets`: projected on tensor and swept axes, isotropic on a generated 2D mesh); used for **the start mesh only**, in ID2's prediction rounds (at most 12, each on the refined mesh, its edges predicted afresh), ID1's guard and AD12's fallback unchanged; every later pass is Unit 17's online layer. (b) Also every later proposal (Unit 17's `prediction` hook). | **(a)**, ruled as recommended: ID8 says the model "predicts the first mesh" and replaces ID2's a priori field; ID9 measures it against 17's online layer, which (b) would remove. A new `PolicyOptions::first_prediction` hook carries it; the existing `prediction` (AC17's test hook) is unchanged. |
| LB5 | Model constants (LQ2's "fixed hyperparameters"). | (a) **Trees:** squared loss, 200 trees, depth 4, shrinkage 0.1, at least 20 samples per leaf, the start the label mean; split candidates at up to 64 cut points per feature (the training values sorted, cut points at equal ranks, duplicates removed); the best split by the largest decrease in squared error, ties to the lower feature index, then the lower cut. Categories as their index. **Ridge (the ID9 baseline):** the continuous features standardized by the training mean and deviation, categories one-hot, the intercept unpenalized, λ = 1 in standardized units, solved by Cholesky. (b) Tuned per training set. | **(a)**, ruled as recommended: fixed, deterministic, single-threaded training; tuning needs a validation split the families are too small for. |
| LB6 | Files and the manifest. | (a) **Adaptation history** `.nithist` and **model** `.nitlearn` in io's container (content kinds 9 and 10; versioned, checksummed, atomic; doubles by bit pattern); the **manifest** a text file (`nitcad-learn-manifest 1`, the schema version, then one line per history: its file digest, FNV-1a 64 of the whole file, and its path), lines sorted by digest, duplicates refused; the training set is the histories' samples in manifest order, within a history in pass then edge order. At most 20 000 samples per history, by a fixed stride (every k-th, k = ⌈n / 20 000⌉), so one large run does not swamp the family. (b) The manifest in the container too. | **(a)**, ruled as recommended: a manifest is written and read by people choosing devices; its digests make it exact. The stride is deterministic, so "the same manifest always gives the same samples" holds. |
| LB7 | The held-out split (ID10: by device). | (a) One history file is one device. `nitcad learn train <manifest> <model> [--holdout <history>...]`: the named histories held out; without `--holdout`, every fourth history in manifest order, starting with the fourth (at least one). Refused (`invalid_input`): no training or no held-out device, a held-out path not in the manifest. (b) A random split with a seed. | **(a)**, ruled as recommended (deterministic; AC19 names its 2 held-out devices). The model is trained on the training devices only and the held-out metrics stored describe that model (it is not retrained on all). |
| LB8 | "Out of the training range" (ID8). | (a) The model stores each continuous feature's training minimum and maximum and the categories seen; an edge is out of range when a continuous feature is outside [min − 5% of the range, max + 5% of it] or a category is unseen (a material name not in the model's list); **a device is out of range when more than 5% of its start mesh's edges are**: the run then uses ID2's a priori field, says so in the report and on the console, and the identity still carries the model (the run asked for it). (b) Per-edge fallback. | **(a)**, ruled as recommended: a per-edge mixture of two predictors would be neither; the device is the unit of validation (ID10). |
| LB9 | The deck. | (a) `[adapt] history = <file>` saves the adaptation history (opt-in, any policy); `predictor = <file>` loads a model (needs `policy = intelligent`, refused otherwise with its line), its digest and versions entering the identity (`adapt.predictor` appended, as AD8's settings); absent: Unit 17's runs bit-identical (AC19). (b) `history` only through the CLI. | **(a)**, ruled as recommended: the deck is the run's whole input. |
| LB10 | ID10's held-out metrics "rank correlation, recall/precision of predicted refinement, passes saved". Passes saved needs adaptive runs on the held-out devices, which training cannot do without their builders. | (a) The file stores the held-out Spearman rank correlation of ℓ̂ against ℓ (trees and ridge), and recall and precision of ℓ̂ ≥ 0.5 against ℓ ≥ 0.5 on the held-out devices' first-pass edges; passes saved is measured by AC19's runs and recorded here, not stored. (b) `nitcad learn train` runs the held-out decks' adaptations. | **(a)**, ruled as recommended; **confirmed by the owner (2026-11-13: "LB10 confirmed, choose recommended")** (it narrows what ID10 says the file stores; no criterion changes: AC19 measures passes saved either way). (b) would make training a simulation campaign. |
| LB11 | Where the history is gathered. | (a) The runner (`deck::run`), through Unit 17's `on_mesh` callback: each solved pass's features computed as it is solved, the labels at the end against the final mesh; the history holds the deck text, the run record (Unit 42 provenance), the schema version, per pass the device and mesh digests, the features and labels. (b) Inside `solve::adapt`. | **(a)**, ruled as recommended: the solve layer stays free of files and deck text (Unit 17's division). |
| LB12 | AC19's family and gate cost. | (a) 8 tensor MOSFETs as deck texts (AR6's tensor MOSFET with L_g, t_ox and N_A varied), `policy = intelligent`, adjoint-weighted, tolerance 1e-2, I_D at V_g 0.5 and 1 V; a manual Release gate `[.ac19]` like `[.ac5]`. The determinism, versioning, identity, fallback and bit-identity clauses are fast cases on small devices. | **(a)**, ruled as recommended. |

**Slices.**
1. Features and labels (LB2, LB3), the history file and its io (LB6), the deck's `history` (LB9, LB11), the manifest
   and `nitcad learn manifest` (LB6); AC19's "the training set built twice from the same manifest is bit-identical".
2. Trees and ridge (LB5), training with the device split (LB7), the `.nitlearn` file with its versions, ranges and
   held-out metrics (ID10, LB8, LB10), `nitcad learn train`; AC19's model bit-identity and version refusal.
3. The predictor in the loop (LB4, LB8), the deck's `predictor` and the identity (LB9), the report; AC19's family gate,
   out-of-range fallback and bit-identity without `predictor`; AC21's predictor case; AC12's ID8 and ID10 mutants;
   ARCHITECTURE 6.12.4 and 11.3.


**Amendment at build (LB9, 2026-11-13).** The history is an output, so it is written by `[output] adapt_history =
<file>` beside `adapt` and `adapt_meshes` (AD8's convention), not by an `[adapt] history` key; `predictor` stays an
`[adapt]` key (it is an input). Nothing else in LB9 changes.

**Implementation record.**

*Slice 1 — features, labels, histories, manifests (LB2, LB3, LB6, LB9, LB11).* `solve/adapt_policy.hpp`:
`feature_points` (ID2's junction, interface and depletion points with their spans and sizes, and the analysis's
largest junction and gate bias differences), split out of `predict_sizes` without a change to it (the fast adapt suite
and `[adapt_intelligent]` unchanged: 28 cases pass; 3 pass and 3 fail as expected). `results/learn.hpp`: schema 1's
17 features, `HistoryPass`, `AdaptHistory`, `Manifest`, `TrainingSet`. `solve/adapt_learn.hpp`: `material_key` (a
category from the material's kind, eps_r and, for a semiconductor, Eg0, Nc300 and Nv300 as exact hexadecimal
floats: materials carry no names, and two decks naming one material differently give one key), `device_digest`,
`edge_features` (nearest junction, interface and contact points by an x-sorted scan, ties to the lower index),
`realized_levels`, `HistoryRecorder`. Feature 10 is ID2's per-point largest junction difference (the diode's points
−1 to 0.6 V give 1 V, not the range 1.6 V). `io/learn_files.hpp`: the `.nithist` file (content kind 9), `file_digest`,
`make_manifest`, the text manifest, `build_training_set` with LB6's stride. `deck::Deck` keeps its text (the history
records it); `[output] adapt_history`; `deck::run` gathers the history through `on_mesh` (after `adapt_meshes`'
writer when both are given). `nitcad learn manifest <manifest> <history>...`. Tests: `[adapt_learn]` (features on
the diode, labels against a refined mesh and against itself, a run's history saved, loaded and saved again
bit-identical, the manifest the same in either input order, the training set built twice equal, the cap, and the
refusals: a duplicate, a digest that does not match, another schema, unsorted lines); the deck case (the history's
deck text, run identity and per-pass mesh digests, the CLI's manifest, its usage and a non-history refused).

**LB13 (found at build, slice 2; ruled as recommended).** ID10 says the model records "the build that trained it
(Unit 42 provenance)", and AC19 that "the trained model" built twice is bit-identical. A Unit 42 provenance holds
the run's start time and machine too, so two trainings would differ. (a) **The model keeps the build's fields only**
(version, revision, dirty, build type, compiler, Eigen version); the run's (start, wall time, OS, CPU, MKL, cache
sizes) are left empty. (b) Keep them and compare models without them. **(a)**, ruled as recommended: the build is
what ID10 names, and the file itself is then reproducible.

*Slice 2 — the model (LB5, LB7, LB10, LB13; ID9, ID10).* `solve/adapt_learn.hpp`: `cut_points`, `fit_boosted_trees`
(histogram split search over LB5's cut points, one byte per bin; the left child grown first; leaves the mean
residual times the shrinkage), `predict_trees`, `fit_ridge` (its 20-odd columns solved by an in-place Cholesky:
Eigen is private to linalg), `rank_correlation` (Spearman, average ranks), `default_held_out`, `train_predictor`
(trees and ridge on the training devices only, both measured on the held-out devices; ranges and categories over
the training samples; the manifest digest over the entries' digests, not their paths), `predict_levels`.
`results::LearnedModel`, `TreeNode`, `HeldOutMetrics`; the `.nitlearn` file (content kind 10) refusing another
format or schema with both versions named, and a tree whose children do not follow their parent or whose feature is
out of range as corrupt. `nitcad learn train <manifest> <model> [--holdout <history>]...` (a held-out path matched
as written or as the same file), printing the metrics and the importances. Tests: the parts exactly (cut points,
rank correlation with ties, the default split, trees learning a step that ridge cannot, the first split, a tie to the
lower feature, `min_leaf`, ridge recovering a linear law), training by device (sides, counts, trees beating ridge on
a step, recall and precision, trained twice equal and its files bit-identical, loaded equal and predicting the
same, another format and another schema refused by name, another content kind, the split's four refusals, a
held-out device's labels not reaching the trees), the diode's histories end to end, and the CLI (training with a
named held-out history, an unknown one, usage).


*Slice 3a — the predictor in the loop (LB4, LB8, LB9; AC19's fast clauses, AC21's predictor case).*
`PolicyOptions::first_prediction`: the start mesh's prediction rounds take it in place of ID2's field (before the test
hook `prediction`, which still drives the later proposals when set: AC21's fixture uses both); every later pass is the
layer's own. `solve::learned_targets` (LB4), `out_of_range` and `out_of_range_fraction` (LB8), `stamp_predictor` (the
identity continued over the model file's digest, format, schema and use; settings `adapt.predictor.digest_high`,
`_low`, `.format`, `.schema`, `.used`). `[adapt] predictor = <file>` is loaded at build, so `nitcad check` refuses a
model of another format or schema with its line, and needs `policy = intelligent`. `deck::run` computes the start
mesh's (the deck's own, AD10) features against the model's materials; above 5% out of range it leaves ID2's field and
records `out_of_range`. `results::AdaptReport` gains the predictor's path, digest, status and fraction (io keys 25–28,
written only with a predictor, so Unit 17's reports are unchanged and still load); `nitcad run` and `check` print
them. Tests: the range rule and the targets exactly (the 5% margin both ways, unseen and absent materials, NaN, the
5% device limit, l^ = 1, 0.4 and the 0.5 threshold); the deck case (`adapt_history` alone leaves the identity and the
deterministic report bit-identical to a run without it; a model used, its digest in the report and its file; two
models two identities; a device of another material 100% out of range, its passes' meshes those of the run without
the predictor, the identity still carrying the model with `used = 0`; the console lines; another format refused at
build); AC21 with a learned model given through `learned_targets`, the certificates bit-identical to basic's.

**AC19 (measured 2026-11-13, Release; `[.ac19]` in tests/deck/cross_run_learning_test.cpp, 18.7 min).** The family
(LB12): tensor MOSFETs from 50 nm × 20 nm uniform starts, source and drain 300 nm, I_D at V_g 0.5 and 1 V, V_d
0.05 V, intelligent, adjoint-weighted, tolerance 1e-2, `max_passes` 10. **ID2's a priori first mesh** (each run's
history saved): m0 (L_g 400 nm, t_ox 4 nm, N_A 1e17) 10 passes, 144 144 DOFs (its stops rejected twice, the basic
loop then grew it to 48 048 nodes, and the verification ran on 198 830: 14.5 min of the gate); m1 (500, 5, 2e17) 5
passes, 21 420; m2 (600, 5, 1e17) 5, 19 026; m3 (700, 6, 5e16) 5, 13 542; m4 (800, 4, 2e17) 4, 12 198; m5 (500, 3,
5e16) 4, 15 867; held out: m6 (700, 5, 1.5e17) 5, 10 476; m7 (600, 4, 7e16) 5, 15 561; every run certified, stopped
converged. **The model** (trained on m0–m5, 94 244 samples after LB6's stride; 21 637 held-out samples): held-out
rank correlation **trees 0.399, ridge 0.161** (the trees beat the ID9 baseline); first-pass recall 0.900, precision
0.722. **Learned first mesh on the held-out devices:** m6 4 passes, **12 000 DOFs** (against 5 and 10 476); m7 4
passes, 10 260 DOFs (against 5 and 15 561); both certified. **AC19's family clause is met on m7 and unmet on m6** (a
pass saved, but 1.15× the DOFs); the gate is `[!shouldfail]` with this measurement. **The other AC19 clauses hold**
(fast cases): the training set and the model built twice bit-identical; another format or schema refused; the
identity changes with the model's digest; an out-of-range device flagged and on ID2's field; without `predictor`,
runs bit-identical. **Ruled under the standing instruction, as Unit 17's unmet criteria were (2026-11-13 ruling after
slice 5b):** the m6 DOF miss is a limitation in ARCHITECTURE.md 11.3, the criterion unchanged. Reading: the learned
field covers the refinement the converged runs made (recall 0.90) and so starts nearer the final mesh (a pass fewer on
both), but at precision 0.72 it also refines where m6's own run did not, and m6's a priori run was already the
family's cheapest.

**AC12 for 17b (hand mutants, 2026-11-13).** C:\Users\disha\mut17\run17b.sh in a separate Release tree (r3): one
exact-string change each; `nitcad`, `nitcad_adapt_test` and `nitcad_deck_test` rebuilt; killed if a build or the
`[adapt_learn]`, `[adapt_policy]` or "deck adapt*" cases fail. **Round 1, 22 mutants, 18 killed.** Survivors: ID8's
stride (2 in place of 1: every capped case had a stride of at least 2), LB5's shrinkage dropped, LB5's `min_leaf` at the
split (the node gate alone held the step test), LB4's predictor on every proposal (no case looked at the later
passes). **Closed by tests:** the uncapped set holds every labelled edge and the capped one exactly LB6's stride; a
one-split tree's leaves are 0.1 of their sides' mean residual; a step at 0.9 cannot be cut with 150 a side; with a
predictor and weighting none every pass after the first is the layer's own (marks). **Round 2: all four killed.**
Every 17b row has a killed mutant and none survives: ID8 (7), LB5 (5), ID10 (5), LB8 (3), LB4 (2), LB9 (1). **AC12 met
for 17b.** (A process note: commit 7359b44 captured round 1's first mutant, applied in the working tree while the
runner ran; 4adba48 restores it, and no commit is made while the runner is active.)

**AC13 for 17b (suites, 2026-11-13, at 5cd7008).** Release: 219/219 passed (868 s). Debug: 201/201 passed (980 s).
The AC19 gate is `[.ac19]` (hidden, Release by hand); the 17b cases in the suites are `[adapt_learn]`, the AC21
predictor case and the deck cases. Unit 17 merged to `main` meanwhile (72377ad, PR #81, owner "merge pr 81"); this
branch holds 17's head (6ffad0b), so its PR to `main` carries the 17b commits only.

### Unit 57 — generic fields and equations, completed (proposed 2026-11-15; approved 2026-11-15; reviewed and amended 2026-11-15; started 2026-11-16; on `assemble/generic-equations`)

**Owner request (2026-11-15):** "What next to implement … as per architecture", then "Plan for finish generic
fields and equations (the Unit 40 follow-up)". This unit does what ARCHITECTURE.md 6.12.1–2 and 11.3 list as "Not
yet" after Unit 40: move the model data out of `DriftDiffusion` and retire it as the wrapper, add the derived field
kind with edge locations, and build `EquilibriumPoisson` from the same terms. The branch starts from `main`
`beb494f`. The number is 57: it was freed when the cuDSS plan was dropped (2026-11-15); 56 is reserved by Unit 55.

**Why now:** 11.1's first priority. Every model since Unit 40 (Units 46–54) was added through terms, but also as
members and data in `DriftDiffusion` (header 1,011 lines, source 2,069, about 53 data members the terms read
through `TermAccess`). Density-gradient confinement (11.3, nanoscale) needs a new unknown per carrier. Unit 30's
open natural-test norm needs the layout of each equation. Both are cheaper on a generic system.

**Standing conditions (as Unit 40):** numerically a pure refactor. F, J and every result are bit-identical
(P-core-pinned traces, Unit 41). There is no change to physics, defaults, the deck, file formats, run identities or
checkpoints. Criteria are fixed before code. Release and Debug suites pass at every slice commit with `/W4 /WX`,
plus hand mutants. Nothing is pushed or merged without authorization.

**Left out of this unit:** dropping the insulator nodes' n and p unknowns (FD5: a layout and checkpoint change);
parallel assembly (FD6: not bit-identical); automatic differentiation (6.12.3); the natural-test norm itself
(Unit 30's follow-up, a consumer of this unit); density gradient and any other new physics; the render layer's
\|E\| and \|J\| (Unit 25b; they read results, not the system).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| GE1 | 0 | **Baselines.** (i) Extend the FG1/FG2 recorder's traced entries to the configurations added after Unit 40: Lombardi, MLDA, Altermatt levels, Schottky contacts, gate tunnelling, energy balance, swept and generated 3D meshes, and an `[adapt]` deck. (ii) **Add trace hooks (`detail::observe`) for what slices 3 and 4 change and the trace does not yet see** (review R3): `EquilibriumPoisson::evaluate` and `residual` (x, F, J), the per-mechanism rates (`RateFields`), the carrier temperatures in K, the gate-tunnelling currents and the energy powers. The hooks only copy values when an observer is set, so no number changes. (iii) **Write GG2's decks that do not exist yet** (R4): the repository has only `decks/finfet_3d*.deck`, `decks/schottky_diode_1d.deck` and the four `examples/`. The new ones go in `tests/decks/` (committed, small, each under a minute in Release). (iv) Extend the FG4 benchmark with an energy-balance and a gate-tunnelling device. Record (i), (ii) and the benchmark on `main`, pinned to P cores. | No assembler change. The recorded files stay in scratch (FD8). |
| GE2 | 1 | **The model data store.** The data fixed at creation moves into an internal `SystemData`, held as `shared_ptr<const SystemData>`. That is the per-node, per-edge and per-contact data, the scaled tables, the node blocks' positions `block_`, the edges' `EdgeTerm::ab/ba`, and the gate tunnels and their positions. The data that changes after creation moves into a value-type `SystemState`. That is the contact biases and metal values, the electrode potentials, gate biases, thermal rises, held carrier temperatures and the minority floor. **It also holds the tunnel paths, their layouts, the base triplets and the Jacobian pattern, because `set_paths` rebuilds all four** (`src/assemble/drift_diffusion.cpp:1824-1841`; review R1). Terms read `(const SystemData&, const SystemState&, const Evaluation&)` and no longer `DriftDiffusion`. The private helpers they call (`edge_fluxes`, `thermal_fluxes`, `level`, `interface_at`, `node_storage`, `path_terms`, `gate_tunnel_fluxes`, `surface_gradient`, …) become functions over the data, in the source file of the model they belong to. `TermAccess` is removed. | A copy shares the data and owns its state, so a relocation in one copy never reaches another. `DriftDiffusion`'s public signatures are unchanged in this slice. |
| GE3 | 2 | **The generic system.** A public `assemble::System`: a composition (terms + `DofMap`) over `SystemData` and `SystemState`, with only the members every equation set has: `create`, `dofs`, `unknowns`, `set_bias`, `make_jacobian`, `evaluate` / `residual` (steady and with a time step), `storage`, the small-signal matrix, `bias_derivative` / `bias_rows`, `update_size` / `apply_update`, `stamp_contacts`, `state_from_potential`. The small-signal path's own interface pass becomes a term, active only in the frequency domain. The model-specific members become free functions in topic headers: `currents.hpp` (terminal, edge, conduction and small-signal currents, resolution, flux terms), `thermal.hpp` (contact heat, electrical power, Joule heat, thermal bias), `tunnelling.hpp` (paths, path states, gate-tunnelling currents), `energy.hpp` (energy powers, hold and release temperatures), `charges.hpp` (contact, gate and trap charges, occupancies). **Queries take `(const System&, x)`. Mutators take `System&`** (`set_thermal_bias`, `hold_temperatures`, `release_temperatures`, `set_paths`, `set_minority_floor`; R6). They reach `SystemData` and `SystemState` through accessors whose parameter is a passkey type only `src/assemble` can construct, so the public header exposes no internal type. Callers in `solve` and the tests move to these (GQ1). The trace labels (`detail::observe` tags) keep their names, so the slice-1 and slice-2 traces compare line for line. | Only `solve` and the tests use `DriftDiffusion`; `deck`, `io` and `app` do not. |
| GE4 | 3 | **Derived fields (FD7).** `Field` gains a kind (unknown, fixed, derived) and the locations edge and interface edge besides node. Shapes stay scalar only: no listed field is a vector, so a vector shape waits for its first user (R9). A derived-field table in `assemble` gives a descriptor and an evaluator `(const System&, x) → values` for each existing diagnostic: band edges and quasi-Fermi levels (node), per-mechanism rates and net rates (node), the electron and hole edge currents and flux terms (edge), trap occupancies (interface edge), and lattice and carrier temperatures in K (node). The solve layer fills `results::BiasPoint`, `RateFields` and `BandDiagram` through the table. | Results structs and file formats are unchanged (GQ4). |
| GE5 | 4 | **Equilibrium from the same terms (FD1).** Equilibrium becomes a composition with a ψ-only `DofMap` and the Poisson, contact, oxide, gate and interface sets. It shares `SystemData`'s scaled tables (both already use `make_scaled_device`). **Shared terms:** only the Poisson flux and the gate term, which compute the same expressions at stride 1. **Equilibrium's own term variants** (R2; `src/assemble/equilibrium_poisson.cpp:140-235`), each keeping its present expressions and order of additions: (i) the charge with carriers slaved to ψ (Boltzmann, Fermi–Dirac, incomplete ionization), assigning the residual and the diagonal; (ii) the insulator rows, which assign 0 to the residual and the diagonal; (iii) the interface edges at Fermi occupancy (`InterfaceEdges::equilibrium`), not the steady-state occupancy; (iv) the contact rows, ohmic at ψ₀ and Schottky at ψ_S, whose diagonal is set to 1 in a final pass after every addition (today's post-pass in `evaluate`). Per entry, the order stays: charge → gate → interfaces → edge flux, then the contact diagonals. `solve::equilibrium`'s signature is unchanged (GQ3). | Removes the separate assembler and its own edge type. The saving is about half of `EquilibriumPoisson`'s 442 lines, since the variants remain. |
| GE6 | — | **Order rule.** As FE3: every entry of F and J receives its contributions in the order and with the operands it has on `main`. `/fp:precise` stays, so moving code between functions or files cannot change a value. The traces check this anyway. | |
| GE7 | — | **Commits.** One branch. Commit 0 is GE1, then one commit per slice, each passing every gate on its own. The owner may stop after any slice; slice 1 is useful alone. Then the docs: ARCHITECTURE.md 2.2, 6.12.1–2 and 11.3 updated to "as built" only after the owner accepts the result. | |

**Criteria (fixed before code):**

| # | Criterion |
|---|---|
| GG1 | **Bit-identical numerics:** on every traced entry (Unit 40's nine plus GE1's), F, every Jacobian value, the time-step and small-signal matrices, `bias_derivative`, every Newton count, step cut, λ, arc step, transient step sequence, current, admittance, field and run identity equal commit 0's bit for bit (P-core-pinned). |
| GG2 | **Byte-identical decks:** the `nitcad run` output tables and every written file (`.nitsweep`, `.nitckpt`, `.nitadapt`) of 12 decks are byte-identical to `main`'s. The provenance record (revision, start time, wall time, R4) is excluded: the files are compared record by record through the `io` reader with that record skipped, and the CSV tables line by line. The decks (GE1 iii): diode, oxide with traps, lumped gate, electrothermal, impact ionization, nonlocal tunnelling, Lombardi + MLDA MOSFET, Schottky diode, gate tunnelling, energy balance, swept 3D FinFET (reduced), and an `[adapt]` deck. |
| GG3 | **FD-Jacobian:** every existing FD gate and per-term check passes unchanged; GE5's slaved-charge term has its own FD check under Boltzmann, Fermi–Dirac and incomplete ionization. |
| GG4 | **No slowdown:** FG4's benchmark with GE1's two devices added; `evaluate` and `residual` medians within 2% of commit 0's, or of its measured spread if larger (FD2). |
| GG5 | **Suites:** full Release and Debug suites pass at every slice commit, with no warnings. |
| GG6 | **Mutants:** at least 8 per slice, each caught by the suites or the traces (the report says which). At least these:
  - a term reading the state's old bias after `set_bias`;
  - **the bias state shared between copies by mistake, so that `set_bias` or `set_paths` on one copy changes another** (R10);
  - the pattern kept in `SystemData`, so that a path relocation is ignored;
  - two terms swapped;
  - the frequency-domain interface term active in steady state;
  - a derived field reading the wrong slot;
  - the equilibrium charge term's Fermi–Dirac branch using Boltzmann;
  - equilibrium missing the interface set;
  - the equilibrium contact diagonal set before the edge flux, not after it;
  - a topic function reading the composition's map of another system. |
| GG7 | **Retired:** `drift_diffusion.hpp` and `TermAccess` no longer exist (GQ1 a). `System`'s public header names no model. A test parses the header's member declarations, not its comments (R11), and requires that none is named for a model: no `thermal`, `heat`, `tunnel`, `path`, `energy`, `temperature`, `schottky` or `gate` in a member name. Adding a term changes only its own file and the composition table. |
| GG8 | **Extensibility (GQ5), two-way coupled (R5):** a test-only equation adds one node field χ to the drift-diffusion composition without editing `System` or any existing term. χ is at every node (`DofMap`'s layout, FD5), held at 0 on insulator and contact nodes. Its row is ∇²χ = ψ − ψ_ref. A second test-only term feeds χ back by adding a source k·χ to the electron and hole continuity rows. The test checks: (i) the `DofMap` and pattern grow by one field; (ii) both terms' partials pass the FD check, including the cross blocks (n, p rows in χ; χ row in ψ); (iii) Newton converges quadratically in its last iterations; (iv) at k = 0 the ψ, n, p fields equal the system without χ to 1e-12; (v) at k ≠ 0 on a 1D bar with a manufactured source, the error falls at second order under two refinements. |
| GG9 | **Layering and scope:** terms, data and state stay internal (`src/assemble`); `System`, `Field`, `DofMap` and the derived-field table are public. No new dependency from `physics`, `mesh` or `linalg`. No change to physics, defaults, deck syntax, file formats or run records. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| GQ1 | The fate of `DriftDiffusion` | (a) Replace it with `assemble::System`, migrate every caller in slice 2, and remove `DriftDiffusion`. (b) Keep `DriftDiffusion` as a deprecated alias for one unit. (c) Keep the name for the generic class. | **(a)**: the callers are `solve` (14 files) and the tests (30), all in-repo; an alias would leave two names for one thing. |
| GQ2 | Where the model-specific queries live | (a) Free functions in topic headers (GE3). (b) Members of `System` (today's shape). | **(a)**: `System` stays generic; a new model adds a header, not members. |
| GQ3 | Equilibrium's public type | (a) `System` created with the equilibrium equation sets; `EquilibriumPoisson` removed, its two users in `solve` migrated. (b) Keep `EquilibriumPoisson` as a thin type over a composition. | **(a)**: one numerical architecture (11.2). Its `carriers` and `charge_neutral_potential` become free functions. |
| GQ4 | Derived fields in results | (a) The table feeds today's result structs; formats unchanged. (b) Results carry a generic list of named fields too (a file-format version bump). | **(a)**: no format change in a refactor; (b) is a later unit, if the app or exports need it. |
| GQ5 | The extensibility test (GG8) | (a) Yes. (b) No. | **(a)**: it is the evidence that density gradient can be added without a second architecture. |
| GQ6 | Slice order | (a) As GE2–GE5. (b) Equilibrium (GE5) before the system (GE3). | **(a)**: GE5 is simpler on the finished `System`. |

**Rulings (owner, 2026-11-15: "choose recommended dont start"):** GQ1–GQ6 all **(a)**, as recommended. GE1–GE7
and GG1–GG9 stand as written. Implementation is not started and waits for the owner's word to begin slice 0.

**Review (2026-11-15; amendments applied on the owner's "do necessary, dont start yet"):** the plan was checked
against the code, and 11 findings were folded into the rows above. No criterion was loosened.
- **R1, GE2:** `set_paths` rebuilds the pattern, base triplets and path layouts, so these belong in `SystemState`,
  not the fixed data.
- **R2, GE5:** equilibrium's contact, insulator, charge and interface rows differ from drift-diffusion's at stride 1,
  so they are kept as variants. The "about 440 lines" claim is corrected.
- **R3, GE1 ii:** the trace did not see the equilibrium assembly, the rates, carrier temperatures, gate-tunnelling
  currents or energy powers, which slices 3 and 4 change.
- **R4, GE1 iii and GG2:** most of GG2's decks did not exist, and the byte comparison had to skip provenance.
- **R5, GG8:** χ was one-way coupled, so it did not show the coupling density gradient needs, and "per semiconductor
  node" contradicted `DofMap`'s layout.
- **R6, GE3:** mutators cannot be free functions over `const System&`. They take `System&` and reach the internals
  through a passkey.
- **R7, start condition (below):** PR #85 changes 10 `solve` files with 29 `DriftDiffusion` references.
- **R8, GE2 notes and GG4:** "copies get cheaper" had no user, since nothing in `src` copies a `DriftDiffusion`. The
  claim and GG4's copy report are removed.
- **R9, GE4:** the vector shape had no field, so it was dropped.
- **R10, GG6:** the "data mutated through a copy" mutant was impossible on `const` data, so it was replaced.
- **R11, GG7:** the header scan reads member names, not comments.

**Start condition (R7):** slice 0 starts after PR #85 (Unit 17c) is merged to `main`, with this branch rebased onto
it. If #85 is still open when the owner says start, the owner is asked whether to rebase onto `mesh/derived-goals`
instead.

**Start (owner, 2026-11-16):** "start Unit 57" while #85 was open (its Release jobs green, Debug running). Asked,
the owner first chose to wait for #85, then asked "cant you do it in Unit 57 without rebase". Checked against
`origin/mesh/derived-goals`: **R7 overstated the overlap.** #85's diff adds or removes no `DriftDiffusion` or
`EquilibriumPoisson` reference; of the files it changes only `src/solve/adapt_adjoint.cpp` (4 lines) uses
`DriftDiffusion`. The rest of the overlap is `DECISIONS.md`, `ARCHITECTURE.md` and `CMakeLists.txt`. So slice 0
starts on this branch as it stands (from `beb494f`, no rebase). When #85 merges, `main` is merged into this branch
and the baselines of what #85 changes (the `[adapt]` deck and the `adapt` traces) are recorded again on the merge
before the next slice is compared against them.

**Slice 0 as built (commit `55f7669`, 2026-11-16; this machine: Core Ultra 5 225F).** No assembler change.
- GE1 ii hooks: `equilibrium.evaluate.{x,residual,jacobian}`, `equilibrium.residual.{x,residual}`,
  `rates.{mechanism,electron,hole}`, `energy_powers`, `gate_tunnelling.{ecb,hvb,evb}`; on the solve side
  `result.{electron,hole}_temperature`, `result.ac_{electron,hole}_temperature`, `result.rates` (each `RateFields`
  vector) and `result.tunnelling_{ecb,hvb,evb}`. Non-empty values in the baseline: 27,050 equilibrium Jacobians,
  37,425 assembler rate vectors, 59,446 result rate vectors, 449 carrier-temperature fields, 221 energy-power records,
  220 assembler and 209 result tunnelling records, so every hook is exercised by the suites.
- GE1 i: `nitcad_adapt_test` links the listener. **Traced entries: 79** (every per-file entry of the assemble,
  solve, deck and adapt executables, the six `solve_*` and two `deck_*` Release entries, `adapt_goal`,
  `adapt_intelligent` and `validation`), each pinned to the P cores (0x3C3) in declaration order through a
  scratch `pin.cmd` wrapper (`start "pin" /b /wait /affinity 3C3`; a bare `""` title is mangled by MSYS). **Two
  full traced runs of the slice-0 binaries (A 13 min, B 14 min) are identical on all 79 entries** (4,437,963
  observed values, 181 MB); the files stay in scratch (FD8).
- GE1 iii: the 12 decks in `tests/decks/` (all converge; the reduced swept FinFET 61–67 s pinned, the others under
  0.7 s). The 1D junction decks use an erfc donor profile on a constant acceptor body: the deck tests'
  `diode_device` (windows `x=-1:1um`, `x=1um:3um`) turned out to give an ohmic, fully compensated bar (its current
  is linear, 7.3 mA/cm² per V), so it was not copied. GG2's comparison is the hidden `[.deck_compare]` case
  (`tests/deck/deck_compare_test.cpp`). **Two runs of the slice-0 `nitcad.exe` agree on all 38 output files.**
- **GG2 and the adaptation report (ruled 2026-11-16, owner "choose recommended": the times are excluded):** the `.nitadapt` report stores measured times (per-pass
  seconds, the cost model fitted to them, the verification and coarse seconds). Between the two identical runs the
  raw report differs in 102 bytes, against about 20 (the provenance) for a sweep file, so GG2's byte comparison
  cannot hold for it with only the provenance skipped. The comparison clears those times as it clears the
  provenance (as ruled; GG2 reads "provenance and, in a `.nitadapt` report, its measured times"); the `[adapt]`
  deck stays in GG2. Slice 1 waits for the owner's word to start.
- GE1 iv: FG4's benchmark gains `Gate tunnelling: assembly cost on the MOS capacitor` (1D and 2D); EC19 is the
  energy-balance device. Baseline, two runs pinned to CPU 0, `evaluate` median in µs (run 1 / run 2): legacy 1D
  diode 12.0/12.1, 2D MOSFET 2642/2494, 3D generated diode 753/761, Zener 12.5/11.7, electrothermal bar 212/213
  (Fermi–Dirac 333/340), energy balance 1D 396/408, 2D 5466/5849, 3D 3489/2771, gate tunnelling 1D 48.1/45.3,
  2D 171/157. The run-to-run spread reaches 26% on EC19's devices, so FG4 is judged as in Unit 40: alternating
  pinned runs, minimum per device.
- Suites: Release 226/227. `adapt_intelligent` failed AC11's timing check once under `-j 6` (a transposed solve
  0.24 ms against one factorization's 0.19 ms on a 73-node diode pass) and passes alone; it is load-sensitive,
  not changed. Debug is left to CI (owner's standing rule).

**Slice 1 (GE2) as built (2026-11-16; owner "start slice 1"; on this branch as it stands, #85 still open).**
- `src/assemble/system_data.hpp` (internal): `SystemData` holds what `DriftDiffusion`'s private section held,
  less the state: the per-node, per-edge and per-contact data, scaled tables, composition, nested types and the
  terms' Jacobian positions. `SystemState` holds the gate nodes (their bias), `psi0`/`n0`/`p0`, the metal values,
  contact biases, thermal rises, held carrier temperatures, the minority floor, and the tunnel paths, layouts,
  base triplets and pattern (R1). `DriftDiffusion` keeps `shared_ptr<const SystemData>` and an owned
  `unique_ptr<SystemState>`; its copy shares the data and copies the state. Its public signatures are
  unchanged; the former inline accessors are out of line.
- Terms take `(const SystemData&, const SystemState&, const Evaluation&)`; `Term`, `Terms`, `Composition` and
  `GateTunnelJob` are `SystemData`'s nested types, and a composition entry's `active` takes both. The helpers
  the terms call are `SystemData` members, defined where they were; those that read the state take it
  (`gate_tunnel_flux(es)`, `thermal_ohmic`, `add_gate_shift`, `assemble_heat`) and the locators take the
  pattern. `TermAccess` is removed; `DriftDiffusion` befriends `SystemData` for `of` / `state_of` alone.
- **As-built note on R1 and the positions.** GE2 lists `block_` and `EdgeTerm::ab/ba` as data, but `set_paths`
  moves every position into the new pattern (`block_`, the edges', the surface stencils', impact ionization's,
  the interfaces', the gate tunnels'). The positions stay in `SystemData`, and `set_paths` moves them in a copy of
  the data (copy on write, `make_shared<SystemData>(*data_)`), so a copy that shares the old data keeps its own
  pattern and positions. `set_paths` runs once per path freeze, so the copy is not on a hot path.
- Tests: two new cases check that a copy owns its state (R10): `band to band assemble: a copy owns its bias and
  paths` (set_bias and set_paths on a copy leave the original's residual, pattern and Jacobian values bit for
  bit; an assigned system takes the copy's) and `electrothermal assemble: a thermal bias moves the sinks of its own
  copy`. Callers outside `src/assemble` are unchanged; the two internal tests rename `TermAccess::` to
  `SystemData::`.
- **GG1:** the 79 traced entries of the slice-1 binaries are byte-identical to slice 0's run A (P-core pinned).
  The trace ran on the build before the long added lines were wrapped (whitespace only) and before the two
  tests were added.
- **GG2:** the 12 decks run with the final `nitcad.exe` agree with slice 0's run A on all 38 files
  (`[.deck_compare]`).
- **GG4:** FG4's benchmark, slice-0 and slice-1 binaries alternating, pinned to CPU 0, minimum per device. Two
  runs each of every `[.assembly_benchmark]` case: every device within 2% of slice 0 on `evaluate` and
  `residual`, or within slice 0's spread, except 1D drift-diffusion + Canali (+3.1%, spread 1.2%), so EC19's case
  was run four more times each: that device then -0.4% / -0.3%; the 3D devices +3.1% to +8.2% on `evaluate`
  against slice 0's spread of 4.2% to 9.9% (in the first two runs 3D drift-diffusion was -0.0%, 623.3 against
  623.2 us); every other device within 2% or faster. Gate tunnelling 1D / 2D: +0.1% / -0.6%. Passed under GG4's
  spread rule.
- **GG6, 9 mutants, all caught** by the `assemble.*` and `solve.*` ctest entries (control 47/47):
  M1 set_bias's changes discarded (terms read the old bias): 31 entries fail; M2 the state shared between copies
  (`shared_ptr`, copy constructor aliases it): `assemble.band_to_band` (the new copy test); M3 `set_paths` does not
  move the positions (the pattern alone relocated): `assemble.band_to_band`, `assemble.terms`, `solve.btbt`,
  `solve.rates`; M4 `set_paths` moves the positions in the shared data (`const_pointer_cast`, no copy):
  `assemble.band_to_band`; M5 SRH and Auger swapped in the table: `assemble.drift_diffusion`, `assemble.terms`,
  `solve.heavy_doping`, `solve.transient`, `solve.small_signal`; M6 electrothermal set_bias leaves the metal value:
  5 entries; M7 set_thermal_bias writes nothing: first survived (no test moved a sink after creation), caught by
  `assemble.electrothermal` after the thermal-bias copy test was added; M8 hold_temperatures keeps no
  temperatures: `assemble.energy_balance`, `solve.energy_balance`; M9 set_minority_floor a no-op:
  `assemble.drift_diffusion`. The slice-2 to slice-4 mutants of GG6's list (frequency-domain term, derived slot,
  equilibrium) belong to those slices.
- **GG5:** full Release suite (`-j 6`, `/W4 /WX`, no warnings) **227/227 pass** in 3,193 s on the final tree.
  Debug is left to CI (owner's standing rule).
- **Incident (slice 0's scratch):** a leftover first version of slice 0's `run_decks.sh` (with `start ""`,
  which MSYS turns into `\\`) was still running from 12:00 and showed "Windows cannot find '\\'" once per deck as
  the owner pressed OK. Stopped at 13:48. It could only rewrite the console logs (`.out`) of run A, which GG2 skips;
  the data files matched.

**Slice 2 (GE3) as built (2026-11-16; owner "start slice 2"; on this branch as it stands, #85 still open).**
- `include/NiTCAD/assemble/system.hpp` (renamed from `drift_diffusion.hpp`, so the header comment keeps its
  history): `assemble::System` with GE3's generic members only, plus the counts and step data every equation
  set has (`node_count`, `contact_count`, `stride`, `time_scale`, `trap_slots`, `storage_width`, `TimeStep`) and
  `bias_nonlinear` (the companion of `bias_derivative`'s state overload). `DriftDiffusion` is gone; every
  reference to it in `solve` (10 files) and the tests (31) now names `System`.
- Topic headers, each a list of free functions over `(const System&, ...)` or, for the mutators, `System&`:
  `currents.hpp` (edge, terminal, conduction and small-signal currents, `CurrentRow`, the resolution, the flux
  terms), `thermal.hpp` (`electrothermal`, `set_thermal_bias`, the sinks, `heat_scale`, contact heat, electrical
  power, Joule heat), `tunnelling.hpp` (`tunnelling`, `trace_paths`, `set_paths`, `paths`, `path_states`,
  `paths_agree`, `gate_tunnelling`, `GateTunnellingCurrents` and their currents), `energy.hpp`
  (`energy_balance`, `EnergyPowers`, `energy_powers`, hold, release, held), `charges.hpp` (trap occupancies,
  contact, gate and interface trap charges). **As-built placements GE3 does not name:** `rates` with
  `RateMechanism` and `MechanismRate` in a sixth header, `rates.hpp`; `band_edges(const System&, x)` in the
  existing `band_edges.hpp` (beside `BandEdges`, with `System` declared, not included); `set_minority_floor` /
  `minority_floor` as free functions in `system.hpp` after the class (a Newton setting, not a model). The two
  internal constants `impact_current_floor` and `matters` are `SystemData`'s (no caller outside `src/assemble`).
- **Passkey (R6):** `System::data`, `state` (const and not) and `replace_data` take a `const SystemKey&`;
  `SystemKey` is declared in the public header and defined in `src/assemble/system_data.hpp` with a private
  constructor that only `SystemData` can call (`SystemData::of`, `state_of`, `replace`). `System` no longer
  befriends `SystemData`. The public header names `SystemData`, `SystemState` and `SystemKey` only as
  declarations.
- **The bodies:** the model queries are `SystemData` members taking the state (`terminal_currents(st, x)`, ...),
  defined where they were; `src/assemble/topics.cpp` holds the free functions, each one line through the
  passkey. Three are free functions in their model's file: `terminal_current_resolution` and `energy_powers`
  call the system's own `residual` (so its trace hooks fire as before), and `set_paths` replaces the data
  (copy on write, GE2's as-built note). The expressions are unchanged: the bodies were moved by rewriting
  `data_->` and `state_->` only.
- **The frequency-domain term:** small_signal_matrix's interface pass is the term
  `Terms::InterfacesFrequency` ("interface traps, frequency domain"), moved with its expressions and order. A
  composition-table row now has a `frequency` flag; the composition keeps such terms in a list of their own
  (`frequency_terms`), so the steady list (`Terms::list`, the per-term checks, `rates`) is unchanged.
  `Term::add_frequency` adds to the complex matrix (nothing by default); the term's steady `add` is a contract
  failure. `Composition::assemble(..., frequency)` leaves out, in the frequency domain, the steady terms of a
  domain a frequency-domain term has (the interface edges), and `assemble_frequency` runs after the storage, as
  the pass did. The electrothermal table has no such term (the model refuses traps and interface charge), so
  its small-signal matrix is as before.
- The relocation loop's stand-in in `tests/solve/btbt_test.cpp` gives the same five functions as free
  functions (found by ADL), as `solve::detail::solve_relocating` now calls them.
- **GG7:** `system: the public header names no model and the old type is retired` (`tests/assemble/fields_test.cpp`):
  `drift_diffusion.hpp` does not exist; no `TermAccess` in the internal headers; the header (comments blanked)
  has no `DriftDiffusion`; in `class System`'s body every declared name (an identifier followed by
  `( ; = { , ) [`: functions, data members, nested types, their members and the parameters) contains none of
  GG7's eight words; and the parse finds the 17 members it must (so an empty parse cannot pass).
- **GG1:** the 79 traced entries of the slice-2 binaries (P-core pinned): 77 byte-identical to slice 1's
  recorded traces. The other two, `assemble.band_to_band` and `assemble.electrothermal`, hold slice 1's two copy
  tests, added after slice 1's trace was recorded; traced again with slice 1's final binaries they are
  byte-identical to slice 2's. The trace ran before the GG7 test was added (it calls no assembler).
- **GG2:** the 12 decks run with the final `nitcad.exe` agree with slice 0's run A on all 38 files
  (`[.deck_compare]`).
- **GG4:** FG4's benchmark, slice-1 and slice-2 binaries alternating, pinned to CPU 0. After two runs each,
  five 1D and 2D devices were 2.7% to 3.9% slower than slice 1's minimum and outside its two-run spread, so
  `[.assembly_benchmark]` of the solve tests ran four more times each (the deck test's gate-tunnelling cases
  were within 2%). Minimum of six per device: every device within 2% of slice 1 or within slice 1's
  six-run spread. The largest are +2.6% (1D energy balance and 2D drift-diffusion + Canali `evaluate`, 2D
  MOSFET `residual`; spreads 8.1% to 10.8%), 2D drift-diffusion `evaluate` +1.9%, and 3D drift-diffusion
  -3.0%. Passed under GG4's spread rule.
- **GG6, 11 mutants, 10 caught** by the `assemble.*` and `solve.*` ctest entries (control 47/47):
  M1 the frequency-domain interface term active in steady state (also in the steady list): 12 entries (its
  steady `add` is a contract failure); M2 the steady interface term kept in the frequency domain:
  `assemble.small_signal`, `solve.small_signal`; M3 the frequency-domain terms never run: the same two; M5 a
  topic function (`terminal_currents`) reading the first system it was called with, not its argument: 24
  entries; M6 `set_minority_floor` writing a copy of the state: `assemble.drift_diffusion`; M7 `set_paths`
  keeping the old data (positions not moved): `assemble.band_to_band`, `assemble.terms`, `solve.btbt`,
  `solve.rates`; M8 `hold_temperatures` holding in a copy: `assemble.energy_balance`, `solve.energy_balance`;
  M9 `set_thermal_bias` writing a copy: `assemble.electrothermal`; M10 `release_temperatures` a no-op: the two
  energy-balance entries; M11 a model-named member (`tunnelling`) added to `System`: `assemble.fields` (GG7).
  **M4, the frequency-domain terms added before the storage instead of after it, survives** the suites and
  the traces (`assemble.small_signal`, `solve.small_signal`, `solve.transient` byte-identical). It is
  equivalent in every traced use: s is purely imaginary there (i omega t0), so the storage adds only to
  imaginary parts, and the two orders compute the same sums (a + b = b + a; the real parts get no storage
  term). It would differ for a real s, which no caller passes. Not counted; 10 caught meets GG6's 8.
- **GG5: not completed.** The final tree builds with no warnings (`/W4 /WX`, Release). The full Release suite was
  started and stopped by the owner after 5 entries (all passed; "stop the test and push"); the branch is pushed
  for CI, whose Release and Debug jobs run the suites. Before the stop, the assemble.* and solve.* entries
  passed 47/47 (GG6's control on the final sources).

**Risks:** the bit-identity surface is the largest since Unit 40. Moving the bias and path state out of the system
could change which values a copy sees (GG6's first three mutants target this). Slice 2's caller migration touches
about 44 files (more if #85 adds callers). Equilibrium's variants (GE5) keep a second set of row expressions, so the
saving there is structural more than in lines. Trace and deck runs need P-core pinning (Unit 41).

---

## New questions the advisor raised (not yet in `ARCHITECTURE.md`)

None of these is decided. Each needs an owner position before the build scaffold fixes the flags.

| # | Question | Needed before |
|---|---|---|
| N1 | Runtime library: pin `CMAKE_MSVC_RUNTIME_LIBRARY` to match `x64-windows-static-md` (dynamic CRT, `/MD`). It is unstated. | Unit 1 |
| N2 | Floating-point model: `/fp:precise` vs `/fp:strict`, and `/arch:AVX2` or not. This affects the tolerance-based reproducibility in 6.8. | Unit 1 |
| N3 | Exception flag: `/EHsc` is needed because Eigen and the standard library can throw. State it. | Unit 1 |
| N4 | Windows headers: define `UNICODE`/`_UNICODE`, use `/utf-8`, and define `NOMINMAX`/`WIN32_LEAN_AND_MEAN` before `windows.h` reaches any code. | Unit 1 (core), L7/L8 (app) |
| N5 | Whether both Debug and Release builds run the tests in CI. | CI unit |
| N6 | Licence. NiTCAD stays MIT (`LICENSE`), as the owner stated with Unit 18 (2026-10-08). Eigen (MPL-2.0) is compiled into the binaries: a release must carry its notice and MPL source terms. Intel MKL (Unit 18) is an optional, separately licensed external dependency (Intel Simplified Software License): NiTCAD loads it at run time from a path the user gives and never links or ships it; bundling MKL with a release would need its licence and notices handled explicitly. | release notices before release |
| N7 | Whether clang-format, `/analyze` or clang-tidy are wanted at all, given "no speculative files". | owner |

## Inconsistencies the advisor found in `ARCHITECTURE.md` (not yet fixed)

These are the advisor's reading; I have not re-checked each against the file. I believe the first three are correct from
my own edits.

1. **Unit 1 gate list.** Section 11 says Unit 1 waits for "V1, V2 and V3". V4 also exists, and V3 (scaling) is consumed at Unit 7, not Unit 1.
2. **V2 row.** V2 was run under `-std:c++latest`, which Q6 rejects. V4 covers the decided flag, but the V2 row reads as if it were current.
3. **Q10 vs section 8.** Q10 says "Unit 1, or when CI is created", while section 8 keeps CI out of the structure until approved. One of them should win.
4. **Warning policy.** Section 9 says warnings are errors, but V4 used `/W4` without `/WX`. The exact flags for Unit 1 are not stated.
5. **Thread count** appears both as a backend option and as a separate item in 6.10 (item 3 and item 5). Minor duplication.
6. **D5 nesting** is listed as provisional in two places (14.1 and a "D5 (part)" row in 14.2). Reference it from one place.

## Suggested review order

1. Q1, Q2, Q4 (small and independent; unblock Unit 1 and Unit 2).
2. Q3, then N1–N4 (these shape the build scaffold).
3. Q5 (before Unit 3).
4. Q10, N5 (when CI is created), Q11 (when 14.52 ships).
5. N6, N7 (when convenient).
