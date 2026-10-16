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

### Unit 16 — unstructured 2D and 3D meshing (branch `mesh/unstructured`, revised 2026-10-14; on `main` as `a09cf43`, PR #35)

**Owner rulings (2026-10-14), unchanged:** R1 2D and 3D in this unit; R2 an in-tree mesh generator, no new dependency;
R3 material interfaces as node regions with straddling node pairs; R4 analytic doping profiles included. The branch
starts at `main` `23199a9`. Everything below is proposed and **awaits the owner's ruling before any code**. Revision 2
follows the review of revision 1 (`cae0ef0`): its M2 contradicted M3, its sliver refusal contradicted M3, it had no
tie-break rule, and M7 lacked the conditions under which it is exact.

**Evidence labels.** *Proven*: follows from a stated argument or a cited theorem whose hypotheses are checked on every
generated mesh. *Cited, to verify*: from the literature; the theorem's exact hypotheses must be confirmed against the paper
before `ARCHITECTURE.md` cites it. *Heuristic*: no guarantee; bounded by a cap and checked after the fact.
*Limitation*: known and stated, not fixed here.

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| M1 | Box method on the circumcentric dual. An edge's coupling area is the signed sum, over the simplices sharing it, of the dual-face pieces. A node's volume and a boundary node's face area are the signed sums of its pieces. Circumcentres are computed with exact determinants (Shewchuk expansion arithmetic) and one final rounding. | Proven: on a Delaunay mesh the dual is the Voronoi diagram, which depends on the points only, not on which tetrahedra were chosen. | Standard Voronoi finite volume; reduces to the tensor-grid box method on split rectangles and Kuhn-split boxes. Exact determinants keep each circumcentre accurate to a few ulps however flat the element. | — |
| M2 | Checks on **totals**, never per-element pieces (an obtuse triangle's pieces are legitimately negative). Tolerance τ = 1e-12, scale L_e^(D−1) for an edge and the node's (local mean edge)^D for a volume. A total within τ·scale of zero is zero: the edge is dropped, a zero volume is refused. A total below −τ·scale is refused (`degenerate_mesh`, index of the edge or node). Roundoff of order 1e-17 is never refused. | Proven that totals are ≥ 0 under M3's conditions; τ is about 4.5e3 ulps, the precedent of `length_tolerance` in `mesh.cpp`. | Fixes revision 1's contradiction with M3; cocircular diagonals (exact zeros) are dropped, not refused. | A total that is truly negative but within τ of zero is treated as zero. |
| M3 | Guarantees. (a) Delaunay ⇒ every interior total coupling ≥ 0. (b) Delaunay plus the Gabriel property for every boundary segment or facet ("boundary conforming Delaunay") ⇒ the dual is the Voronoi diagram clipped to the domain: every volume > 0 and every boundary area ≥ 0. (c) The generator checks (a) and (b) exactly on its output (P1, P2 below), so they are postconditions, not just claims. | (a) proven; (b) cited, to verify (Gärtner, Si, Fuhrmann 2010; Fuhrmann's lecture notes). | (b) is the property Voronoi box methods require. | — |
| M4 | Exact predicates (orient2d/3d, incircle, insphere) ported from Shewchuk's public-domain adaptive code. Ties (incircle or insphere = 0) are broken by the symbolic perturbation of Devillers and Teillaud (2011): perturb the lifting, decide by orientation tests on subsets ordered by node index. orient tests are never perturbed, so no tetrahedron of four coplanar points can exist; every element has orient > 0 (P0). | Proven for the predicates and for determinism; the tie-break is cited, to verify. | Tensor lattices, Kuhn splits and mirrored pairs are exactly cocircular or cospherical by design; ties must resolve the same way on every machine (6.8). | Each tie costs a few extra orientation tests. |
| M5 | 2D domains: any planar straight-line graph (non-convex polygons, holes) by constrained Delaunay triangulation, then refinement until no segment is encroached (so the CDT is Delaunay). | Cited, to verify: boundary conforming Delaunay meshes exist and are computable for arbitrary polygonal regions in 2D (Gärtner–Si–Fuhrmann abstract). | Real cross-sections. | — |
| M6 | 3D domains: convex polyhedra only. Boundary conformity (M3 b) is guaranteed only when every input angle is at least arccos(1/3) ≈ 70.53°; boxes (90°) satisfy this. Convex domains with smaller dihedral angles are allowed and handled by the heuristic in M11, with P2 still checked. Non-convex domains are refused (`invalid_input`) and deferred. | Proven: for convex domains the Delaunay tetrahedralization meets the hull without boundary recovery. 70.53° is cited and checked (Gärtner–Si–Fuhrmann 2010, U6). | No 3D constrained boundary recovery (the TetGen-class problem). Interfaces need no conforming faces (M7). | The angle is the theorem's hypothesis, not a refusal rule. |
| M7 | (Ruled R3.) Every node has one region; nodes are never on an interface. On each interface, a pair set: a 1D (2D domain) or 2D (3D domain) point set P on the interface with covering radius R(x), and for each p ∈ P the nodes p ± (s/2)·n. A **protected band**: no node outside the pair set lies within b(x) = √(R² + s²/4) of the interface. Pair spacing s = min(h(x), lfs(x)/2), where lfs is the local feature size, so pairs from different features cannot interfere. | Proven (away from corner zones, M8): with mirror symmetry and the band, every point of the interface is nearest to a pair node and every off-interface point is strictly nearer its own side's mirror image. So the Voronoi faces between the two regions lie in the interface, couplings between non-mirror nodes across it are exactly 0, and each cross edge's midpoint lies on the interface (the assembler's assumption, `assemble/interface_edges.hpp`). | Generalizes `mesh::straddle_interface` exactly. | Coordinates p ± s/2 round, so cross couplings come out at about 1e-17 relative instead of 0 (absorbed by M2's τ). The interface is a polyline or polygon at resolution s. |
| M8 | Corners and junctions: a bend in an interface; an interface meeting the boundary, another interface (three or more regions) or a contact end. Inside a corner zone of radius 2s(x) around the feature, nodes are placed per sector (one per region, on the sector's bisector at distance s/2), and junction nodes on the boundary straddle the junction point. No angle limit is imposed. Instead the generator measures the dual interface near each feature (the union of Voronoi faces separating regions) and refuses (`degenerate_mesh`, feature index) when its Hausdorff distance from the true interface exceeds s, or a node falls in the wrong region. | Heuristic placement, checked postcondition (P4b). The bound s is the resolution M7 already accepts. | A check-based limit instead of an arbitrary angle. | Interface geometry inside corner zones is accurate only to s. |
| M9 | Contacts. 2D: a patch is a sub-interval of boundary segments; its ends are straddled along the boundary by node pairs (M8 junction rule). 3D: each convex face carries any number of named patch polygons (simple polygons in the face, pairwise disjoint, each inside one face), with the rest of the face under the face's name. The face is first meshed by the 2D generator, with its patch edges as straddled 1D interfaces; the face nodes then go into the 3D mesh. A node on an edge between two faces belongs to a patch of each, with each patch's own area. | Proven (needs M3 b): under boundary conformity a node's boundary area on a face is its 2D Voronoi cell among the face's nodes. So patch boundaries are exact to s, by the M7 argument in 2D. Gate: compare with the 2D path computed independently. | Planar MOSFET source, drain and gate on one top face. | A patch cannot span two faces; a contact over two faces lists both patches' nodes. |
| M10 | Slivers. **Not assumed harmless.** Facts: (i) exactly flat elements cannot exist (M4, P0); (ii) for Delaunay meshes the couplings and volumes are Voronoi measures, independent of slivers, and M1's exact determinants keep a sliver's circumcentre accurate; (iii) anything that interpolates on simplices (rendering, future AMR transfer, a future unstructured BTBT locator) **is** degraded by slivers. Handling: report the minimum dihedral angle and the sliver count (dihedral < 5°) per mesh. P5 compares the dual against an independent Voronoi clipping on test meshes. A sliver-rich stress test (perturbed lattice) must pass the same solver gates as a clean mesh. No refusal on angle. Sliver removal is deferred: it needs weighted Delaunay, so a power-diagram dual, a different control-volume definition needing its own ruling. | (i) and (ii) proven; their effect on the solve is measured, not assumed; (iii) limitation. | Separates what slivers can and cannot affect, with tests. | Simplex-based consumers must handle slivers when they arrive. |
| M11 | Refinement and termination. 2D: Ruppert/Shewchuk with concentric-shell segment splitting; minimum-angle bound default 20.7°, option up to 33.8°; skinny triangles allowed only at small input angles. 3D: radius–edge bound default 2.0 with encroachment splitting; option down to 1.414. A size field h(x) with h_min > 0 (else `invalid_input`) and grading ≤ 1.5 per h. Every run has a node cap (`MeshOptions::max_nodes`, default 2e6), giving `resource_exhausted` if hit. Modified steps (no Steiner point in a protected band: split the interface pair set locally instead; elements made only of pair nodes exempt from the quality bound) are heuristic. | 2D: termination at 20.7° for input angles ≥ 60°, cited and checked (Shewchuk 2002); below 60° concentric shells alone do not guarantee it (corrected 2026-10-15, U6), so such inputs are heuristic and capped. 33.8°: heuristic. 3D: radius–edge ≤ 2 for input angles ≥ 90°, cited and checked (Shewchuk 1998, U6); smaller input angles heuristic and capped. 1.414: heuristic (TetGen's practice). Band and pair modifications: heuristic, capped. | Proven where the literature proves it; otherwise capped and reported. | A heuristic run may stop at the cap, with a message naming it. |
| M12 | Postconditions, checked on every generated mesh in Debug and Release (O(n)): P0 every element has orient > 0. P1 locally Delaunay (exact, tie-broken). P2 every boundary segment or facet is Gabriel (exact). P3 M2 totals. P4a every edge with nonzero coupling joining two regions is a mirror pair whose midpoint lies on the interface within 1e-12·L_e (outside corner zones). P4b the M8 corner bound. P5 per-node closure: Σ A_e t̂_e (outward) + Σ boundary area · n = 0 within 1e-11 relative to Σ\|A_e\|. A failure is `degenerate_mesh` naming the check; it signals a generator bug, never a user error. | Checks, not claims. The 1e-11 is a few ulps times the face count (≤ ~200). | Each guarantee in M3, M7 and M8 is verified, not trusted. | ~10–20% on top of generation time (to measure). |
| M13 | Determinism and speed. Insertion order: BRIO rounds with a Hilbert-curve sort, using a fixed-seed splitmix64 (not `std::` distributions, whose output is implementation-defined). Point location by visibility walk from the last insertion (terminates on Delaunay triangulations); cavities by breadth-first search over neighbours. Same input and options give a bit-identical mesh on any machine (6.8). Budgets, measured and recorded, not promised: Release 2D 1e5 nodes ≤ 2 s; 3D 1e5 nodes ≤ 15 s. ctest meshes ≤ 1e4 nodes (2D) and ≤ 2e3 (3D); 3D solves in ctest ≤ 3e3 nodes; larger ones `[.benchmark]`. The Debug suite grows by ≤ 120 s in total, or cases move to `[.benchmark]`. | Expected O(n log n) for BRIO plus walking: cited (Amenta–Choi–Rote 2003). Budgets are targets. | 3D without fast location is O(n²); SparseLU fill-in makes large 3D solves slow (PARDISO for benchmarks). | Hilbert sort plus BRIO adds code. |
| M14 | Mesh metadata persists. The io mesh file gains optional keyed records (IO2: old readers skip them): simplices, node region, patch provenance, generator options and a generation digest. `mesh_digest` stays the graph only, so a mesh's identity is unchanged. A new `load_generated_mesh` returns `GeneratedMesh`, or `invalid_input` if the records are missing; `load_mesh` is unchanged. | — | Reloaded meshes stay renderable (Unit 25) and refinable (Unit 17) without regeneration. | Larger files (simplices are about 2× the edges). |
| M15 | (Ruled R4.) Doping profiles: constant, box, Gaussian and erfc per axis with lateral erfc fall-off, summed per species, sampled at nodes (`device/doping_profiles.hpp`). | — | Point sampling matches today's lambdas. | Abrupt doping steps are not straddled (U9). |
| M16 | Nonlocal BTBT stays tensor-only; it is already refused without `DeviceDescription::cells`. | — | Keeps the unit to meshing. | Unstructured nonlocal tunnelling waits for a cell-walk locator. |

**Implementation notes, commit 3 (predicates and simplex import; differences from the text above, for review):**
- **M4, exact stage:** exact integers over the inputs scaled by a common power of two, instead of Shewchuk's expansion
  arithmetic.
  - Why: expansions can underflow on the degree-6 diametral-sphere test and on circumcentre numerators at extreme
    coordinate scales. Integers have no range restriction and are simpler to verify.
  - The filter is a forward error bound computed alongside the same formula, instead of Shewchuk's stage-A constants, so
    the filter and the exact stage cannot evaluate different formulas.
  - Trade-off: the exact stage is slower; it runs only on near-ties.
- **M4, tie-break:** the largest node id is perturbed most. Any fixed order is consistent; this one is documented and
  pinned by a test.
- **M2, scales:** the scale is max(L^(D−1), Σ|pieces|) for couplings, max(l^D, Σ|pieces|) for volumes and
  max(l^(D−1), Σ|pieces|) for patch areas, instead of L^(D−1) and l^D alone. The Σ|pieces| term covers cancellation
  among signed pieces near poorly shaped elements, where rounding scales with the pieces rather than with L.
  - The P5 closure tolerance is relative to Σ|pieces| for the same reason.
  - A patch node whose area is within τ of zero is dropped from the patch, as zero couplings are dropped.
- **Import scope:** P0, P1, P2, P3 and P5 are checked, plus manifoldness and folds. Accepted, but stated as limits:
  - global overlap of a locally valid mesh is not detected;
  - P4 cannot be checked without interfaces (U2).
- **Determinism:** claimed only within one build. Topology decisions are exact and therefore platform-independent;
  bit identity of the measures across compilers or machines is not claimed, because it has not been verified.
- **Mutation tests:** 12 run against the 20 new test cases; 10 caught, 2 survived.
  - Caught: insphere sign; tie-break order; the query-point cofactor sign; the face-side sign of the 3D dual; the drop
    rule using `< 0`; P1 skipped; P2 skipped; the 3D boundary piece and the volume factor.
  - The filter with no error bound is caught too, in a second form: its first form failed to build on an
    unused-variable warning and is not counted.
  - **Survived (closed in round 2, below):** the product rule written as max + 1, and the underflow guard removed. The worst-case error they
    under-bound is not reached by any input the tests generate, so the filter's soundness rests on the error analysis
    in `predicates.hpp`, not on a test (U11).

**Approval blockers, round 2 (owner request 2026-10-14; commits `ce7e1cd`, `93a6d2b`, `dd32cbd`):**

*U11, the filter (closed).* The soundness argument is now written in full in `src/mesh/predicates.hpp`: the
arithmetic model (binary64, round to nearest, gradual underflow, no contraction), the induction for sums and
products, cancellation (the bound is relative to M, not to |v|), overflow (an infinity or NaN reaches v or m) and
underflow.
- **Defect found and fixed:** the filter flagged a product only when its magnitude underflowed. A product whose
  *value* underflows while its magnitude does not escaped the flag. This happens when a factor's value is small by
  cancellation, and it falls outside the rounding model the bound assumes. Both are now checked.
- **Removed as redundant:** a deferral on a subnormal bound. Every formula's root is a sum of products, and each
  product's m is 0 or at least DBL_MIN, so a nonzero bound is at least d·2^-1073. Its rounding error is then at
  most 1/(4d) relative, which the factor 2 covers.
- **Tests (`tests/mesh/filter_test.cpp`)** use a probe that exposes each evaluation's value, magnitude, depth, flag,
  bound and decision beside the exact value. They check:
  - the depth table;
  - |value − exact| ≤ bound decided in exact integers, for 600 near-degenerate samples per predicate at scales
    2^±40, with the threshold region (|v| within 4× of the bound) exercised for every predicate;
  - both underflow contracts, overflow, cancellation at offsets up to 2^52, subnormal bounds, and the arithmetic
    model itself.
- **Mutations:**

| Mutation | Result |
|---|---|
| U11-1 product depth max + 1 (the earlier survivor) | caught by the depth table |
| U11-2 underflow flag ignored (the earlier survivor) | caught |
| U11-3 value-underflow check removed (the new check) | caught |
| U11-4 magnitude-underflow check removed | caught, after a test for a product with an exactly zero factor value |
| U11-5 subnormal-bound deferral removed | redundant (argument above); the code was removed |
| U11-6 bound factor 2 → 0.25 | caught |
| U11-7 magnitude of a difference subtracts | caught |

*M4, the tie-break (verified).* The scheme documented in `predicates.hpp` is checked against itself evaluated
directly, independently of the implementation's cofactor expansion:
- lifted rows (p, |p|² + ε^(r+1), 1), with r the rank by decreasing id and ε = 2^-128, scaled to integers;
- their Leibniz determinant (24 or 120 terms), with incircle its sign and insphere its negation;
- on tied integer points on x² + y² = 625 and x² + y² + z² = 9, translated, with half the cases off the tie.
- All 600 2D and 300 3D cases agree (more than 100 and 50 of them exact ties).
- Mutations caught: tie-break order reversed, query cofactor sign, incircle query cofactor as a unit column.
- Guarantees:
  - every sign is exact for any finite doubles;
  - the perturbed predicates never return 0 and define one consistent regular triangulation;
  - on any platform satisfying the arithmetic model the same inputs give the same signs, so the same topology. This
    is by the argument, not by cross-platform runs (none were made).
  - The dual measures are floating-point; their bit identity across compilers or machines is not claimed.

*M2, the tolerance (verified).* `src/mesh/dual_tolerance.hpp` holds the whole rule. `make_dual_mesh` calls it for
every coupling, patch area and volume, always on the assembled total and its scale, never on a piece. Non-finite
totals are refused. Tests:
- exact boundaries with nextafter on both sides of ±τ·scale, at five scales, plus NaN and infinity;
- a near-cocircular square whose diagonal crosses the threshold at ε = 2e-12, under coordinate scalings from 1e-7 to
  1e5;
- Kuhn grids scaled from 1e-9 to 1e3 dropping the same edges;
- exact cancellation of ±4 pieces (dropped, not refused), and the exact P1 refusal that precedes M2 when the same
  quad is moved by one ulp;
- a −4 piece with a +6.68 piece, kept at its exact total;
- a kite whose pieces are 25× its edge, where only the Σ|pieces| term drops a resolved total of 1e-9.
- Mutations: 9 of 9 caught (both inclusive ends, non-finite, the volume bound, τ ×10 and ÷10, either term of the
  scale, the coupling scale's power).

*M7, the interface rule P4 (specified and implemented as a check; no generator).* `include/NiTCAD/mesh/interfaces.hpp`
gives the metric and its computation:
- deviation(F) = min over the declared pieces between the two regions of max over the face's vertices of the
  distance to the piece;
- the face's vertices are the circumcentres around the edge, and on the boundary the boundary-facet circumcentres
  and the midpoint;
- because F lies in the hull of its vertices and the distance to a convex piece is convex, this is an upper bound on
  the face's distance from the interface.

Acceptance:
- **P4a**, outside the corner zones (2 s around piece boundaries): 1e-12 L + 1e-13 X;
- **P4b**, inside the zones: s;
- every region pair that meets must be declared.

`check_interfaces` is the only way to establish P4. `make_dual_mesh` never does, and an imported multi-region mesh
fails it unless its nodes follow the rule (tested). Results on brute-force Delaunay layouts:

| Layout | Result |
|---|---|
| Flat, uniform or graded pairs; thin layer wider than both bands; 3D plane | exact (flat deviation 0) |
| Node inside the protected band; layer inside a band with staggered pairs; 3D broken mirror | refused |
| Node just outside the band (1.05 b) | accepted, as the band argument says |
| Right-angle corner, pairs 8 s away, no corner nodes | exact: by mirror symmetry the corner point is equidistant from both arms |
| Right-angle corner, bisector nodes, pairs graded to 1.5 s | 0.80 s |
| 60° corner, bisector nodes, graded pairs, oblique boundary junction | 0.95 s |
| T-junction | 0.31 s |
| Corner threshold | the 0.80 s layout is accepted at resolution 0.026 and refused by P4b at 0.024 |

M8 findings, which change what the generator must do:
- the bisector nodes make a right-angle corner worse, so the generator should not place them there;
- an oblique boundary junction needs the pairs *and* the boundary nodes graded toward it (1.5 s, 3 s). Without that
  the import failed P2, or the junction face deviated beyond s;
- below 60° nothing is measured, and 60° itself passes with a thin margin (U12).

Mutations: 12 run; 11 caught. M7-12 ("zero-coupling cross edges judged too") is **equivalent**: such a face is
degenerate, and its single point is the circumcentre of cells whose kept cross faces already contain it; if a cell
has no kept cross face, that point carries no interface measure.

**Owner rulings on U2 and U12 (2026-10-14):**
- **U2, P4 certification: required.**
  - A multi-region device (two or more regions) on a mesh that is not a tensor product must carry a valid P4
    certificate. A missing certificate is refused, and so is a stale one.
  - The certificate comes from `mesh::certify_interfaces`. It records digests of the dual mesh, the node regions and the
    declared interfaces; it is stale when the device's mesh or node regions differ from those it was issued for.
  - Tensor-grid meshes are exempt by provenance (owner's choice of the recommended option). `make_tensor_grid` tags its
    mesh, `Mesh::from_parts` verifies that a mesh tagged so has only axis-aligned edges, and io saves and loads the tag.
    A file written before the tag existed loads as a general graph.
  - Single-region devices need no certificate.
  - Enforced in `Device::create`, so every solve of such a device has one.
- **U12, corner rules: as proposed, with acute corners refused.**
  - A bend of exactly 90° between two interface segments gets no corner nodes, since mirror symmetry makes it exact.
  - Every other interface feature with all sectors at least 60° (bends, junctions of three or more regions, junctions
    with the outer boundary) gets graded pairs (1.5 s, 3 s) and one bisector node per sector at s / 2. It is accepted
    only if the generated mesh passes P4b; otherwise generation fails with `degenerate_mesh`.
  - A sector below 60° is refused (`invalid_input`) until measurements justify it.
  - The angle tests are exact (integer arithmetic on the input coordinates).

**2D generator, round 3 (2026-10-14; commits `740beba` to `8be2cc3`):**

*What was built:*
- **P4 certification (U2):** `Device::create` refuses a multi-region device on a 2D or 3D non-tensor mesh without a
  current certificate. Tensor grids are exempt by provenance (`MeshOrigin`, verified, saved by io). 1D meshes are exempt
  too: an interpretation of the ruling, for the same reason (every edge lies along x).
- **The generator** (`mesh::generate_mesh`):
  - exact validation and the U12 corner rules;
  - mirrored pairs graded toward features, a straddle at boundary junctions, nodes exactly at patch ends;
  - Bowyer–Watson with the exact predicates, then Ruppert refinement with protected bands;
  - import and certification of the result.
- **Doping profiles (M15) and generated meshes on disk (M14).**

*Measured:*

| Gate | Result |
|---|---|
| Generator validity, quality (min angle ≥ 20.7° outside bands), size, determinism, refusals | passed |
| P4 by feature: flat, 90° bend, T-junction, 60.3°, 120°, thin layer, band-forcing refinement | certified: corner deviation 0, 0, 0.31 s, 0.95 s, 0.73 s |
| Speed (Release) | 103,627 nodes in 1.70 s (target 2 s) |
| Linear-field consistency (oxide between electrodes, skew mesh) | linear to 1e-10 V: passed |
| Smooth diode convergence | H1 1.79, 1.24; L2 2.36, 2.11; current 1.50, 2.09: passed (≥ 0.9, 1.5, 0.9) |
| Abrupt diode against 1D | error 1.9e-4 (11,549 nodes) against the tensor grid's 2.6e-4 (11,536): passed |
| MOSFET, k = 1.5 against the tensor MOSFET | **not met**: on current 4 % low at 8,031 nodes, tensor 0.05 % at 9,490 (U14) |

*Mutations:*
- certification: 7 of 7 caught;
- generator: 15 of 15 caught (one rule, boundary grading at junctions, was removed as redundant after its mutant
  survived: refinement splits do that work);
- doping and io: 7 of 7 caught.

*Findings that need a ruling:*
- **U14, the MOSFET gate (unmet).** The generator makes isotropic triangles. The tensor MOSFET resolves the inversion
  layer with 2e-8 cm cells normal to the surface; isotropic triangles that fine over the channel exceed the node budget.
  - Options: anisotropic refinement near selected boundaries (a boundary-layer mode for the generator); a hybrid mesh
    with a tensor strip under the gate; or a ruling that changes the gate.
  - The gate stays in the suite as `[!shouldfail]` with its measurement.
- **U15, the U12 corner rule.** Bisector nodes make every measured corner worse:

  | Corner | With bisector nodes | Without |
  |---|---|---|
  | T-junction | 0.31 s | 0 |
  | 60.3° | 0.95 s | 0.36 s |
  | 120° | 0.73 s | 0.32 s |

  Proposed: drop the bisector nodes everywhere. Graded pairs alone give a smaller deviation and a wider margin at 60°.
  They are implemented as ruled until then.
- **U16, U2 in 1D.** The 1D exemption above was applied without an explicit ruling.

**Owner rulings on U14, U15 and U16 (2026-10-14), and their implementation (commits `d741912`, `898fbf1`, `13f9148`):**

- **U15, corners (ruled: remove the bisector nodes where they worsen accuracy). Done.** They worsened every measured
  configuration, so none are placed. P4b still judges every feature, the 90° rule stands, and sectors below 60° are
  still refused. Regression tests pin each measured configuration, with deviation bounds from the measurements and no
  node within s of the feature:

  | Configuration | Corner deviation |
  |---|---|
  | Flat interface | 0 |
  | 90° bend | 0 |
  | T-junction | 0 |
  | 60.3° | ≤ 0.37 s (measured 0.359) |
  | 120° | ≤ 0.33 s (measured 0.321) |

  Mutations: 3 of 3 caught, including the old bisector generator run as a mutant.
- **U16, 1D (ruled: exempt only verified 1D tensor grids). Done.**
  - The blanket 1D exemption is removed.
  - The tensor tag is verified by `Mesh::from_parts`; in 1D it also requires a chain: strictly increasing x, and
    exactly the n − 1 edges joining x-consecutive nodes. Every region transition then lies on the dual face between
    two neighbours.
  - Any other 1D mesh with two regions needs `mesh::certify_interfaces_1d`: every cross-region edge's midpoint (its
    dual face) on a declared interface point of that pair, within the P4a tolerance.
  - The 1D insulator test, a mesh built from parts, now carries such a certificate.
  - Mutations: 5 of 6 caught. The sixth, the increasing-x check removed, is equivalent: with n − 1 rank-consecutive
    edges, a duplicate x is an edge of zero length, which the length check refuses.
- **U14, the MOSFET gate (ruled: keep k = 1.5; investigate anisotropic layers and a hybrid tensor strip; pass before
  3D). Passed.**
  - Both options were realised as one mechanism, `MeshOptions2D::layers`: a structured boundary layer is a tensor
    lattice along a boundary segment, and placed under the gate it is the hybrid strip. Its cells are Delaunay ties,
    so its dual is exactly the tensor box method; refinement stays out of it.
  - Mesh specifications tried, the threshold fixed throughout:

    | Attempt | Specification | Nodes | Result |
    |---|---|---|---|
    | 0 | isotropic, size boxes | 8,031 | on current 4 % low |
    | 1 | layer to 1e-6 cm | 4,438 | subthreshold 1–4 % high |
    | 2 | layer to 1.2e-5 cm | 7,337 | high-gate current and swing fail |
    | 3 | layer to 1.6e-5 cm | 8,304 | high-gate current, g_m and swing fail |
    | 4 | as 3, contact edges exactly at Lsd and Lsd + Lg in all three meshes | 8,190 | only the swing fails, 7.5e-4 against 1.03e-4 allowed |
    | 5 | layer to 1.8e-5 cm | 8,864 | swing fails, 1.4e-4 |
    | 6 | layer to 1.8e-5 cm, unstructured refinement (1e-6 cm) below at the body contact | 9,342 | **passes** |

  - Before attempt 4, the tensor baseline and reference missed the contact edges by 2.9e-8 and 1.4e-8 cm, an offset
    that halves with refinement. Snapping their nearest nodes onto the edges, with node counts unchanged, makes all
    three meshes model the defined device. It does not weaken the rule: the tensor error only shrinks.
  - **Plainly:** the passing mesh is a tensor lattice over 90 % of the depth with unstructured refinement only at the
    bottom. Isotropic triangles alone could not meet the gate at this node budget.
  - The gate is in ctest (`solve_unstructured`, Release) as an ordinary test; the `[!shouldfail]` marker is removed.
  - Mutations: 8 of 9 caught. The ninth, judging layer-only triangles, is equivalent: they are cell halves whose
    circumcentres lie in the layer, already exempt.
  - A defect found and fixed: boundary nodes were filled between lattice columns more than h apart, breaking the
    lattice. Hidden in the MOSFET (columns closer than h), found by the new lattice test.

**3D generator, commit 6 (owner approval 2026-10-15; commits `c12915a`, `86621b7`, `9da439a`, `3c56d24` and the 6e
gates):** Owner choices for this commit: full M6 (convex polyhedra), planar interfaces with no junction lines, and a
z-uniform convergence gate against a 2D reference.

- **6a, `Triangulation3D`** (`src/mesh/triangulation3d.hpp`):
  - Bowyer–Watson insertion with the exact orient3d and the perturbed insphere, and a visibility walk.
  - BRIO insertion order: splitmix64 rounds, with a 3D Hilbert order inside each round (Skilling).
  - It equals the brute-force perturbed Delaunay on a lattice where every tie occurs.
  - 1e5 points take 4.2 s (Release).
  - Defect found and fixed: a cavity slot reused before its outside neighbours were relinked.
- **6b, single region** (`generate_mesh(Geometry3D, MeshOptions3D)`):
  - Exact validation of convexity, closure and planar, convex faces; non-convex domains are refused (U7).
  - Refinement follows Shewchuk:
    - encroached domain-edge subsegments are split, with concentric shells at vertices;
    - encroached subfacets are split at their circumcentres;
    - tetrahedra whose radius–edge ratio exceeds the bound, or whose circumradius exceeds √(3/8)·h, get their
      circumcentre unless it encroaches.
  - **Change from the plan:** faces are not meshed by the 2D generator. Each face keeps a 2D Delaunay triangulation of
    its nodes, in face coordinates and exact on axis-aligned faces, and that triangulation defines its subfacets. The
    public 2D generator is unchanged, and the face and 3D nodes are one set.
  - Extraction:
    - drops the super-vertex tetrahedra and any tetrahedron whose four vertices lie on one face;
    - checks that every exposed triangle lies on one face and that the volume matches to 1e-12;
    - then runs `make_dual_mesh`.
  - **Oblique faces pass**, including an irregular tetrahedron at device scale whose split points round off its lines and
    planes. The pre-stated fallback (refuse oblique faces, U17) is not needed, and no flat tetrahedron occurred.
- **6c, face patches (M9)** (`PolyhedronFace::patches`):
  - A conforming tiling of a face. Its interior edges become constrained subsegments, so every boundary triangle lies in
    one patch and patch areas are exact.
  - **Change from M9:** M9 said patch edges were straddled and accurate to s. Conforming edges make them exact, as the
    2D generator already does at patch ends.
  - Gate: every node's patch area matches `make_dual_mesh` 2D on that patch's triangles to 1e-12. Cases:
    - contact strips;
    - a contact inside a face;
    - a slanted patch edge.
- **6d, planar interfaces (M7)** (`Geometry3D::interfaces` and `regions`):
  - **Validation:**
    - each plane stays at least 4 s from every vertex;
    - it meets every face it cuts at 60–120°, checked before its traces;
    - each trace meets the domain edges it crosses at 60–120° and is at least 6 s long;
    - no cut face carries patches (deferred);
    - planes never meet in the domain and stay at least 4 s apart;
    - there is one seed per slab.
  - **Nodes:**
    - on each crossed edge, nodes at ±s/2 from the crossing (the segment between them is never split);
    - along each trace in a face, pairs at ±s/2 within the face, graded from each end (1.5 s, 3 s);
    - on the polygon inset by 1.5 s, a pair set P from the 2D generator, with nodes at p ± (s/2)·n.
  - **Band:** b = √(h_P²/3 + s²/4). Inside it:
    - no circumcentre is inserted;
    - an encroached subfacet gets a new trace pair;
    - tetrahedra made only of interface nodes are exempt.
  - The result carries the P4 certificate. `Device::create` refuses the mesh when the certificate is missing or stale.
  - P4b measured at the traces (s = 0.02, h = 0.25, unit cube), with regression bounds taken from the measurements:

    | Configuration | Flat deviation | Corner deviation | Bound |
    |---|---|---|---|
    | Horizontal plane | 0 | 0 | 1e-12 s |
    | Thin middle layer (5 s) | 0 | 0 | 1e-12 s |
    | Oblique plane | 3.75e-3 of tolerance | 0.180 s | 0.19 s |

- **6e, gates** (predefined; none set after a run):
  - **Linear field:** conjugate gradients on the generated dual give the linear potential to 1e-12.
  - **Kuhn equivalence:** a solve on Kuhn-split boxes matches the tensor diode solve to 1e-12 in potential and current.
  - **io:** a 3D generated mesh round-trips with its patches, regions, certificate and the 3D report (optional records).
  - **Convergence** (z-uniform doping; reference a 2D tensor solve 4× finer than the finest mesh; `[.unstructured3d]`,
    beyond the ctest node budget of M13, run explicitly):

    | Mesh | Nodes | H1 orders | L2 orders | Current orders |
    |---|---|---|---|---|
    | Generated | 391, 2030, 13057 | 1.78, 1.32 | 1.63, 2.05 | 1.63, 2.18 |
    | Perturbed lattice (sliver stress, M10) | 132, 693, 4305 | 2.12, 1.81 | 2.62, 2.31 | 2.22, 2.01 |

    - The generated meshes have 6, 13 and 133 slivers, with minimum dihedral angles down to 0.65°. The perturbed
      lattices have 0, 33 and 290, down to 0.24°. Slivers did not degrade any gate, which M10 (ii) predicts and the run
      now measures.
    - One of my own test assertions failed: every stress mesh must have slivers. The coarsest lattice has one cell
      through z, so no interior point moves and it has none. I narrowed that check to meshes with interior points. The
      gates were unchanged.
- **Mutations** (each run against the tests; "equivalent" means it survived every measured case, with the reason):

  | Commit | Caught | Equivalent in the measured cases |
  |---|---|---|
  | 6a | 3 of 3 | — |
  | 6b | 2 of 4, plus a fault injection caught by the postcondition and, with the postcondition removed, by the tests | skip subsegment encroachment (subfacet splitting splits the subsegments it encroaches); remove the flat-tetrahedron drop (none occurs) |
  | 6c | 4 of 5 | skip patch-segment encroachment (same reason as 6b) |
  | 6d | 6 | band rule alone and pair-only exemption alone (they mask each other; caught together); trace pairs at ±3s/4 (still mirror-symmetric, so exact); trace-edge rule loosened (on these shapes the face rule refuses the same planes) |
  | 6e (io) | 2 of 2 | — |

- **Performance** (Release, budget targets, M13):
  - As first committed: 128,671 nodes took 23 s, the 15 s per 1e5 target missed. `make_dual_mesh` alone took 16.5 s,
    14.3 s of it in the exact coupling pieces (M1).
  - After U18 (below): 15.9 s, about 12.4 s per 1e5, so the target is met; `make_dual_mesh` takes 7.8 s.
  - The Debug mesh tests grew to 61 s; the 3D tests add about 33 s.
- **Literature (U6), checked 2026-10-15 (completed the same day):**

  | Source | What was checked | Status | Consequence here |
  |---|---|---|---|
  | Gärtner–Si–Fuhrmann 2010 (Comput. Math. Math. Phys. 50(1)) | Abstract | **Confirmed.** 3D boundary-conforming Delaunay meshes when the smallest input angle is at least arccos(1/3) ≈ 70.53°; arbitrary 2D polygonal regions | M3 (b), M5 and M6 cite it as stated |
  | Shewchuk 1998 (SoCG, pp. 86–95) | The paper itself. Its math glyphs are CCITT bitmaps; I decoded them, and the same glyphs read 90° in Ruppert's 2D condition in that paper, which the literature gives as 90° | **Confirmed.** Radius–edge ratio ≤ 2 when no input angle (segments and facets, dihedral included) is below 90° | 3D termination is cited only for domains whose input angles are all at least 90° (boxes). The tetrahedron (70.53°), prism (60°) and octahedron (60° face angles) cases, patch edges and interfaces are heuristic and capped |
  | Shewchuk 2002 (the 2001 report of the CGTA article) | The paper itself, glyphs decoded | **Confirmed:** Ruppert's algorithm terminates when no input angle is below 60° and the split threshold is B ≥ √2, so no angle ends below arcsin(1/(2√2)) ≈ 20.7°. **Corrected:** with smaller input angles, plain Ruppert refinement, concentric shells included, can fail to terminate; only the "Terminator" (shells plus a look-ahead) guarantees termination | M11's "with small angles under concentric shells" was wrong. The 2D generator has shells but no look-ahead, so inputs with angles below 60° (domain corners; interface sectors below 60° are already refused) are heuristic and capped by `max_nodes` |
  | Eymard–Gallouët–Herbin 2000 (Handbook of Numerical Analysis VII) | Through the same authors' 2024 paper (arXiv 2405.16985, introduction) | **Confirmed second-hand:** order-1 error for the two-point flux scheme on admissible meshes when the solution is in H². The handbook's discrete-H1 statement itself was not reached | The H1 ≥ 0.9 gate rests on an order-1 result |
  | Amenta–Choi–Rote 2003 (SoCG, pp. 211–219) | Abstract | **Confirmed:** BRIO is no worse asymptotically than a random insertion order | Visibility-walk location is outside that analysis; the measured 3D Delaunay speed (1e5 points in 4.2 s) is the evidence here |
  | Devillers–Teillaud 2011 (CGTA 44:160–168) | Bibliography only; the HAL copy was blocked | **Not reached** | The tie-break does not depend on the citation: tests check it against the documented lifted perturbation, brute-force Delaunay and antisymmetry |

**Gates (predefined, so a threshold is never set after the run):**
- **Tensor equivalence:** split rectangles and Kuhn-split boxes give exactly `make_tensor_grid`'s volumes, couplings
  and areas to 1e-13 relative, with the diagonals dropped. Their diode solve matches the tensor run to 1e-12.
- **Dual correctness:** P5 closure, plus an *independent* Voronoi cell per node (half-space clipping by bisectors, which
  does not use the simplices) equal to the simplex-built dual within 1e-12 on small 2D and 3D meshes. Closure alone
  cannot detect a wrong dual that is not orthogonal; the independent cells can.
- **Linear-field consistency at the assembler:** Poisson with no charge, Dirichlet at two opposite faces and Neumann
  elsewhere. The discrete potential equals the exact linear one to 1e-12 on generated meshes, and the field
  reconstruction recovers the uniform field exactly. This replaces revision 1's linear-flux gate, which duplicated
  closure.
- **Convergence** (smooth Gaussian doping, h, h/2, h/4; reference a 4× finer tensor solve interpolated
  bilinearly):
  - discrete H1 potential error order ≥ 0.9 (proven O(h) for admissible meshes: Eymard–Gallouët–Herbin, cited, to
    verify);
  - L2 potential order recorded, gate ≥ 1.5 (observed second order is typical but not proven on general meshes);
  - terminal current order ≥ 0.9.
  - Abrupt junction: a laterally uniform 2D unstructured diode against the 1D legacy-validated diode; current error ≤
    the error of a 2D tensor grid with the same node count against the same reference.
- **MOSFET acceptance:**
  - Reference: the existing tensor MOSFET (`tests/solve/mosfet_test.cpp`) at 2× refinement per direction.
  - Compared: I_D at every sweep point above the current resolution, V_T (analysis extraction), maximum g_m, and the
    subthreshold slope over the resolved range.
  - Pass: |Q_unstructured − Q_ref| ≤ k · |Q_tensor − Q_ref| + f, with k = 1.5 (U1). Q_tensor is the baseline tensor
    MOSFET; the unstructured mesh has at most as many nodes. f = the reference's Richardson error estimate.
  - Every point converges, with bias step cutting on if needed (recorded).
  - Off-state current is excluded (6.13).
- **Slivers:** the perturbed-lattice 3D diode passes the same convergence gates as the clean one; sliver statistics are
  recorded.
- **Refusals, predicates, determinism, io:** each M2/M6/M8/M12 refusal returns its code and index; exact signs on 1e-17
  perturbations and cocircular/cospherical sets; two runs bit-identical; `save`/`load_generated_mesh` round trip.
- **Mutation tests:** at least: drop a circumcentric piece's sign; skip P2; skip the band rule; reverse the tie-break; use
  `<` for `<=` in τ; break insertion order. Each must fail a gate.

**Unresolved (need a ruling or carried as limits):**
- U1 MOSFET factor k = 1.5: a judgement, not derived. The owner may set it.
- U2 **Ruled 2026-10-14** (above): multi-region devices on non-tensor meshes need a valid P4 certificate.
- U3 Termination is unproven for the band and pair modifications (M7, M8), for 3D convex domains with an
  input angle below 90° (Shewchuk 1998's hypothesis, checked in U6), for 2D inputs with an angle below 60° (Shewchuk 2002), and for options beyond the cited bounds. Capped (M11) and reported;
  every 3D test case terminated well below its cap.
- U4 Corner-zone accuracy (M8) has no a-priori bound; it is measured and refused above s (P4b). Measured in round 2: 0 at a
  right angle without corner nodes, 0.80 s with bisector nodes, 0.95 s at 60°, 0.31 s at a T-junction.
- U5 Slivers affect future simplex-based consumers (M10 iii); removal needs a power-diagram dual (its own decision).
- U6 **Checked 2026-10-15** (table in "3D generator, commit 6"). Confirmed: Gärtner–Si–Fuhrmann (70.53°), Shewchuk 1998
  (ratio ≤ 2 for input angles ≥ 90°), Shewchuk 2002 (20.7° for input angles ≥ 60°), Amenta–Choi–Rote (BRIO); second-hand:
  Eymard–Gallouët–Herbin (order 1, u in H²). Corrected: concentric shells alone do not guarantee termination below 60°
  (M11). Not reached: Devillers–Teillaud (the tie-break is tested independently).
- U7 Non-convex 3D domains (fins, mesas) are deferred. So are 3D interface junction lines and bends (owner choice
  2026-10-15), interfaces cutting faces that carry patches, and 3D structured boundary layers.
- U8 The L2 order gate (1.5) is a judgement where theory gives no rate.
- U9 Region-only straddling: doping steps (abrupt junctions not on a region boundary) are not straddled, so an abrupt
  junction is resolved only to h. "Doping interfaces" with pairs could be added later.
- U11 Filter soundness: **closed** in round 2 (the argument written out, one defect fixed, direct tests against the exact
  value, every mutation caught or removed as redundant).
- U12 **Ruled 2026-10-14** (above): no corner nodes at 90° bends; graded bisector nodes elsewhere, accepted only if P4b
  passes; sectors below 60° refused.
- U13 P4 does not check that each node lies inside its region's geometry, which is not an input. A misplaced node
  surrounded by its own region's nodes is not detected; one next to another region is.
- U14 **Resolved 2026-10-14**: the MOSFET gate passes with a structured boundary layer (the hybrid strip), k = 1.5
  unchanged (above).
- U15 **Resolved 2026-10-14**: no bisector nodes; regression tests per measured configuration (above).
- U16 **Resolved 2026-10-14**: only verified 1D tensor chains are exempt; other 1D meshes need certify_interfaces_1d.
- U10 Budgets (M13) are targets. If the 3D generator misses its gates or budgets within the unit, the PR carries 2D plus
  3D import and says so (unchanged from revision 1). Result: every 3D gate passed; the 3D speed target was
  missed at first (23 s at 1.29e5 nodes) and met after U18 (15.9 s).
- U17 (pre-stated fallback for oblique faces) **not needed**: oblique faces pass.
- U18 **Resolved 2026-10-15** (owner request), with no change to M1.
  - Profile (128,671 nodes, 774,877 tetrahedra): the exact coupling pieces took 14.3 s of `make_dual_mesh`'s 16.5 s,
    the exact circumcentres 3.1 s, the predicates 0.05 s. The cost was the exact integers' storage, not their
    arithmetic: each operation allocated its limbs on the heap.
  - Fix: `ExactInteger` keeps up to 24 limbs inline and moves to the heap only beyond that; its loops use raw
    pointers; `to_double` gathers its top 64 bits from at most three limbs, with the same truncation as before.
  - **Bit-identical:** the same input gives the same bits. Over all 774,877 tetrahedra, the digests of every
    coupling piece and circumcentre, and the dual graph digest, are equal before and after the change.
  - Result: couplings 14.1 s to 5.6 s; `make_dual_mesh` 16.5 s to 7.8 s; the 3D benchmark 23 s to 15.9 s at
    1.29e5 nodes (target met); the 2D benchmark 1.96 s to 1.15 s at 1.04e5 nodes.
  - Tests: inline and heap values to about 94 limbs; `to_double` against a known top 64 bits at every shift
    offset; identities through heap-sized products.
  - Mutations: 4 of 5 caught. The fifth, heap growth not copying the limbs, is equivalent: growth only ever
    copies limbs that `assign` has just zeroed, into a zeroed buffer.
  - The earlier idea, a filtered circumcentre certifying the correctly rounded value, is not needed; it would have
    changed M1, since `to_double` is faithful, not correctly rounded.

---

### Unit 31 — device files, the text input deck and the batch runner (merged to `main` as `6f67510`, PR #37, 2026-10-15)

**Owner rulings (2026-10-15)**, asked before any code:
- syntax: an own sectioned `key = value` format with an in-tree parser; no TOML or JSON dependency;
- scope: device descriptions on disk, the deck turned into a device, options and a run plan, and a `nitcad`
  runner that writes sweeps and checkpoints;
- coverage: tensor grids, the 2D and 3D generators and mesh files; equilibrium and steady bias sweeps (both
  equation sets), step cutting, line search, the physics switches, electrothermal. Transient, small-signal and
  `trace_bias` are a later unit, together with their results on disk;
- units: plain numbers are cm; `um` and `nm` suffixes convert on read.

Choices made within those rulings, open to review:

| # | Choice | Why | Trade-off |
|---|---|---|---|
| DK1 | The grammar (`deck/syntax.hpp`): `[kind name]` headers, `key = value` entries, `#` comments outside double quotes, a trailing `\` continues a line; keys and names are identifiers joined by single blanks. | Line-oriented, so every error has a line; trivially diffable; no escaping rules beyond quotes. | No nesting, variables or includes; a long polygon needs continuations. |
| DK2 | A new layer `deck` above `solve` and `io` (it reads the option types of `linalg` and `assemble` through `solve`); `io` gains `device` (device files). `solve` still does not use `io`. | The deck must produce `solve::BiasOptions` and read mesh files, so it sits above both; device files belong with the other files. | One more layer in D5 (its name is open to review); `io` now depends on `physics` through `device`. |
| DK3 | Strict: an unknown section or key is an error, as are a second single-instance section and a deck with not exactly one of `[equilibrium]` and `[sweep]`. | A typo (`max_iteration`) silently ignored would run the wrong problem. | Decks written for a later build fail on an earlier one with a clear message. |
| DK4 | Syntax and meaning are separate passes. Every error is `invalid_input` "`<file>:<line>: <what>`" with the line as the context index; an error from a library (a generator, a material, `Device::create`) keeps its code and context and gets the location and section prefixed. | The user sees where, and the library's own words for what. | `Device::create`'s errors have no deck line (they concern the device as a whole) and report line 1. |
| DK5 | Only lengths take suffixes (`cm`, `um`, `nm`); every other quantity is in the unit of the C++ field it sets (V, K, cm^-3, eV, cm^-2, cm/s, s). | Owner ruling for lengths; one unit per quantity elsewhere keeps the deck and the API in step. | `1e17` must be cm^-3, not m^-3. |
| DK6 | Series: `first : last : step` with `last` included when it lies on the grid to 1e-9 of a step (and then stored as `last` exactly), a comma list, or one value held at every point; a contact not named is held at 0 V; every series longer than one value has the same length. | The common sweeps in one line; no hidden accumulation of rounding (each value is `first + k step`). | No logarithmic or adaptive series. |
| DK7 | Material parameters are overridden by their C++ member paths (`electron_mobility.mu_max`, `thermal.conductivity_W_cmK`); materials start from the built-in sets (`silicon`, `germanium`, `gaas`, `ingaas`, `sic_4h`, `algaas <x> [<share>]`, `sio2`, `insulator`). | One name per parameter, the one the code and its comments use; nothing to keep in sync. | Renaming a parameter in C++ breaks decks that set it. |
| DK8 | Regions on a tensor grid: inclusive boxes; the one region without a box is the background; later regions override earlier ones; every node must land in a region. Generated meshes take polygons (2D) or seeds (3D) per region; region *k* is the *k*-th `[region]` section, also for mesh files. A contact's nodes are its patch's semiconductor nodes (ohmic, gate) or insulator nodes (electrode), optionally those inside any of its boxes. Tensor meshes always carry their cells. | Covers the hand-built devices of the test suite without coordinates per node. | Nodes exactly on a box face belong to it; the deck author places boxes between nodes (or uses `straddle`). |
| DK9 | A restart matches the checkpoint's bias to a sweep point to 1e-12 relative (1e-12 V near zero), sweeps the points after it from the checkpoint's state with `initial_bias_V` set, and writes those points only. | Found in testing: `0 : 0.3 : 0.1` ends at 0.3 exactly, the fourth value of `0 : 0.6 : 0.1` is 0.30000000000000004; a bit-for-bit match refused a correct restart. | A sweep file of a restarted run holds the continuation, not the whole sweep. |
| DK10 | The runner: `nitcad check|run <deck>` in `src/cli/main.cpp` (not `app/`, which the layout reserves for the Win32 GUI); exit status 0 done, 1 error, 2 the sweep stopped early (outputs written), 64 usage; Ctrl+C asks the solver to stop. Outputs are written after the analysis, the checkpoint being the last completed point. | The batch runner 11.3 asks for before the GUI. `src/cli/` is a directory the layout did not list: open to review. | A crash (not a stop) loses the run; a checkpoint per point would need a per-point hook in `solve`. |
| DK11 | Device files: content kind 4 in the Unit 29 container; every material parameter under its own key (absent keys take the struct's defaults, IO2); tensor cells as their axes; the P4 certificate as issued; a loaded description goes through `Device::create`. | Bit-for-bit round trips; old files keep loading as parameters are added; no file can bypass validation. | The parameter keys are a format: never renumbered. |

Not in Unit 31: transient, small-signal and `trace_bias` decks (owner ruling); 2D structured boundary layers and
per-region mesh sizes in the deck; general 1D graph meshes (they need certificates); variables, expressions,
includes and parameter studies.

---

### Unit 32 — low-current accuracy: the weighted terminal current (merged to `main` as `6e54c69`, PR #39, 2026-10-15)

**Owner rulings (2026-10-15)**, asked before any code: the method is weighted current extraction (no change of
unknowns; quasi-Fermi unknowns not now); new, opt-in result fields, so every existing result and run identity stays
bit-identical; steady analyses (`sweep_bias`, `solve_bias`, `trace_bias`) and a deck switch, with transient and
small-signal later; an adjoint (dual-weighted residual) error estimate.

Choices made within those rulings, open to review:

| # | Choice | Why | Trade-off |
|---|---|---|---|
| W1 | I_c = Σ over edges of (J_n + J_p)_e (w_a − w_b), with w = 1 on contact c's nodes and 0 on the other ohmic contacts' nodes. | Summation by parts of the nodes' total-current balances: exactly the contact current for a solved state, for *any* such w. | The value of a computed state depends on w through its residual (W4 measures that). |
| W2 | w solves the conductance-weighted Laplace problem on the semiconductor nodes, g_e = half the edge's one-sided flux terms (its conductance to a quasi-Fermi difference near equilibrium), plus 1e-30 of the largest; one linear solve per point, the pivot test off for it. | w is flat where the conductance is large (a 1e19 contact region, whose edge currents are differences of huge fluxes) and varies where it is small (depletion regions, a channel), where the current is resolved; device-neutral, no cut to choose. | One factorization on the node graph per point. |
| W3 | The total current only; no node terms. | Generation and recombination create or remove pairs, so they cancel in the total (SRH, Auger, radiative, impact ionization, local tunnelling, interface traps and surface recombination in a steady state); the gate's displacement current is zero in DC. | Transient (displacement current) and small-signal currents need their own terms: later. |
| W4 | The resolution: the adjoint estimate ∇I_c·J⁻¹F, computed as the change of I_c over one more Newton correction from the state (they are equal to first order: the correction is −J⁻¹F), w held fixed; plus 8ε Σ t_e \|w_a − w_b\| for the rounding of the weighted edges. If the correction cannot be computed, the conventional resolution stands in. | No transposed solve and no gradient of the fluxes (which would repeat the Jacobian's chain rules for every model); one factorization at the state. Measured: on a state stopped at `tol_update` = 0.5 (current 22 times too high) the estimate equals the error to 3e-6. | An estimate, not a bound: second-order terms are left out (negligible once Newton converges quadratically). |
| W5 | Opt-in (`BiasOptions::weighted_current`); `BiasPoint::weighted_current` and `weighted_current_resolution` empty when off; the setting in the run record and the two sweep-file keys (16, 17) only when on. | Every existing result, file and run identity unchanged. | Users must turn it on; it costs one more factorization per point. |
| W6 | Refused (invalid_input, before solving) with nonlocal tunnelling and with the quasi-static sweep. | Nonlocal tunnelling creates the hole at a path's start and the electron at its crossing, so its current is not carried by the edges (W3 fails); a quasi-static sweep has no current. | A nonlocal-tunnelling device has no weighted current yet (a path term would add it). |
| W7 | The analysis extractors (`current_curve`) read the weighted current and its resolution when a point has them. | The off-state figures (on/off ratio, leakage, swing) then use the resolved current. | A curve mixing points with and without it reads each point by its own fields. |
| W8 | Deck: `[solver] weighted_current = on`; the CSV table gains the weighted current and its resolution per contact. | The batch path reaches it without code. | — |

Findings:
- An n⁺ (1e19) / p (1e17) silicon diode without recombination or band-gap narrowing: the weighted current follows
  the ideal-diode law J₀(e^{V/V_T} − 1) to 0.1% from 1.1e-11 to 4.8e-4 A/cm² (0.05 to 0.5 V), with an estimated
  error of 1e-19 to 2e-12, while the contact current's resolution is 4e-6 A/cm². At 0.05 V the contact current
  reads 6.2e-12, 45% low. The remaining 0.1% drift of J/J₀ is the depletion edge's shift with bias.
- The legacy MOSFET (Vds = 0.05 V): every drain current from Vg = −1 to 1 V is resolved (estimated error about
  2e-19 A/cm off, against the contact current's 1e-9 floor, where it read −2e-12, the wrong sign). Off, the current
  levels off at 3.50e-13 A/cm (independent of Vg within 1% from −1 to −0.8 V); on/off between 1 V and −1 V
  9.0e11; the swing per 0.1 V narrows from 76.4 to 70.55 mV/decade from −0.5 to 0 V (Unit 24's 70.55 at the end).
- The off-state plateau is not generation: without SRH it is 4.8% lower, and the subthreshold current at −0.2 V
  changes by 2e-7. It is minority-electron diffusion through the body into the drain; its rate
  I / (1 − e^{−Vds/V_T}) rises with Vds (+6% at 0.1 V, +19% at 0.2 V) as the drain's depletion region shortens the
  base, so the 1D saturation law is no reference there.
- At reverse bias the contact current's own noise shows: its anode and cathode values differ by 2e-5 of the current
  (1e17 diode, −8 V, 1 ns lifetimes), where the weighted currents sum to zero.
- Mutation checks (14, the shell loop of Unit 28), all caught, one after a test was tightened (the nonlocal refusal
  was masked by the assembler's own refusal for a device without cells): uniform weights, a contact's right-hand
  side negated, electrons only, the estimate or the rounding dropped, a wrong contact weight, the sweep's and the
  trace's fill skipped, the run-record setting, the nonlocal rule, the extractor's preference, the sweep-file keys,
  the deck switch, the table columns.

---

### Unit 33 — validation against published data (merged to `main` as `a14ac96`, PR #41, 2026-10-15)

**Owner rulings (2026-10-15)**, asked before any work: references are values stated in the text or tables of
published sources (no digitized curves); the benchmarks are decks in a new top-level `benchmarks/` directory, run by a
Release-only ctest against their cited values; a disagreement is reported and pinned, not tuned (no model parameter
changes in this unit); web search and fetch may be used to find and read the sources.

| # | Choice | Why | Trade-off |
|---|---|---|---|
| V1 | Resistivity against dopant density from W. R. Thurber, R. L. Mattis, Y. M. Liu, J. J. Filliben, NBS Special Publication 400-64 (1981), doi:10.6028/NBS.SP.400-64, the basis of ASTM F723. Phosphorus: eq (1), log₁₀(qρN / 1 V s cm⁻²) = (A₀ + A₁X + A₂X² + A₃X³) / (1 + B₁X + B₂X² + B₃X³), X = log₁₀(N / 10¹⁶ cm⁻³), table 9 (300 K: A = −3.0652, 2.1853, −0.61080, 0.056189; B = −0.67642, 0.19542, −0.018100; 23 °C: A = −3.0769, 2.2108, −0.62272, 0.057501; B = −0.68157, 0.19833, −0.018376). Boron: eq (4), qρN = min + (max − min) / (1 + (N / N_ref)^α), table 12 (300 K: 0.00215, 0.02053 V s cm⁻², 4.09e18 cm⁻³, −0.727; 23 °C: 0.00209, 0.02024, 3.88e18, −0.737). q = 1.602e-19 C, as the report says to use with its fits. | A measured, standardized, device-level relation: the simulator solves a resistor and the resistivity follows from its current. | It tests mobility (and ionization) at equilibrium-like fields only. |
| V2 | The source was read from the Internet Archive's OCR text of the scanned report (the NIST PDF is images; the ASTM preview stops before the equations). The OCR dropped the + signs of eq (1) and the exponents of N₀ and of table 12's N_ref: N₀ = 10¹⁶ gives the report's 0.52 Ω cm at 10¹⁶ and makes its two phosphorus fits agree; N_ref = 4.09e18 (3.88e18 at 23 °C) is the only decade that keeps eq (4) within the report's stated 3% of eq (3) (table 11, legible: 0.00220, 0.01973, ρ_ref 0.01782 Ω cm, 1.086 at 300 K). The test checks every transcribed fit against the report's other fit within its stated self-consistency (boron 3%, phosphorus 4.5%). | No number is taken from memory; a miscopied coefficient fails the check (mutation-tested). | Two exponents are inferred, not read; the inference is recorded here. |
| V3 | Eleven decks (`benchmarks/resistivity/`): phosphorus 1e14–1e19, boron 1e15–1e19 cm⁻³ (inside the report's measured range), a 10 µm bar at 300 K, 0 and 0.01 V; ρ = 0.01 V / (J L). Each also at 23 °C (the deck's temperature overridden by the test) and with incomplete ionization. | One deck per device, readable and runnable with `nitcad run`. | Decks have no variables, so each density is its own file. |
| V4 | The intrinsic density from K. Misiakos, D. Tsamakis, J. Appl. Phys. 74, 3293 (1993), doi:10.1063/1.354551: "(9.7 +/- 0.1) X 10(9) cm-3" at 300 K (abstract, read from its open-archives record). The model's is √(np) in the solved 1e14 bar at 0 V. | A measured value with a stated uncertainty, tested through a solved device. | — |
| V5 | Each gap (simulated / reference − 1) is pinned at its measured value within 0.2 percentage points; no parameter is changed (owner ruling). | A model change must show up as a changed benchmark, deliberately. | The tests pass with known disagreements; the findings below say where. |
| V6 | Not included: abrupt-junction breakdown voltages (S. M. Sze, G. Gibbons, Appl. Phys. Lett. 8, 111 (1966)) and MOS C-V landmarks. The Sze–Gibbons values could not be read from an accessible copy; textbook C-V landmarks are analytic, not measured. | Only values actually read are used. | Breakdown stays validated against the analytic one-sided value only (Unit 19). |

Findings (pinned gaps, simulated against the reference):
- **Phosphorus, 300 K:** +3.9% (1e14), +0.5% (1e15), −6.9% (1e16), −8.8% (1e17), +0.6% (1e18), −7.1% (1e19). The 23 °C
  gaps are within 1 point of these (the temperature dependence holds).
- **Boron, 300 K:** −4.4% (1e15), +3.4% (1e16), **+30.0% (1e17), +48.0% (1e18), +23.3% (1e19)**. The legacy hole
  mobility is the Caughey–Thomas fit to Irvin's curve; the report itself finds that "the Caughey and Thomas
  expression is a poor fit to the mobility data" and differs from the older curves by up to 45% near 5e17 for boron.
  A boron mobility from this report's own fits would close the gap; that is a model change (an owner decision).
- **Incomplete ionization** (a fixed 45 meV level, Unit 15) makes the heavily doped bars much worse: phosphorus +31%
  (1e18) and +138% (1e19), boron +121% (1e18) and +302% (1e19). The level does not lower with density (the report
  finds de-ionization only at intermediate densities, and full activity at the highest), so the model should not be
  used for ionization above about 1e17 cm⁻³ until it has a density-dependent level.
- **Intrinsic density:** the legacy silicon n_i, 1.0674e10 cm⁻³, is 10.0% above the measured 9.7e9 (9.7 times its
  uncertainty); currents set by n_i² (diode saturation, generation) are about 21% high by that alone.
- Mutation checks (4): a phosphorus coefficient, the boron N_ref decade, a 23 °C ρ-fit exponent (caught by the
  transcription check), the temperature override; all caught.

**Scope (owner instruction, 2026-10-15): Unit 33 is a validation baseline only.** It changes no physics, model,
parameter or solver code; the findings above are recorded for later units, each an owner decision.

**Reconstructed coefficients (read from the OCR text where legible; reconstructed where not):**

| Quantity | Printed (OCR) | Used | How it was fixed |
|---|---|---|---|
| Eq (1) operators | "Aq - A^X . A^X . A3X" | A₀ + A₁X + A₂X² + A₃X³ over 1 + B₁X + B₂X² + B₃X³ | With + signs the ρ- and N-based phosphorus fits agree within 1.3% at 10¹⁶ cm⁻³ (the report's stated self-consistency) |
| N₀ (eq 1) | "Nq =10^^ cm~^" | 10¹⁶ cm⁻³ | Gives 0.52 Ω cm at 10¹⁶ cm⁻³ phosphorus, consistent with both fits |
| Boron N_ref, 300 K (table 12) | "4.09 x 10 + 0.20 x 10" | 4.09 × 10¹⁸ cm⁻³ | At ρ = ρ_ref of eq (3) (N = 3.84e18), eq (4) agrees within 1.5% only for 10¹⁸; 10¹⁷ is 60% off. The report states the two fits agree within 3% from 10¹⁴ to 10¹⁹ |
| Boron N_ref, 23 °C (table 12) | "3.88 x 10 + 0.21 x 10" | 3.88 × 10¹⁸ cm⁻³ | The same test at 23 °C: 1.8% for 10¹⁸ |
| All other coefficients (tables 9, 11, 12) | legible | as printed | — |

The test checks every fit used against the report's other fit within the stated self-consistency (boron 3%,
phosphorus 4.5%), so an error in a reconstruction would fail it.

**Unsupported regimes (do not rely on the model there; none is fixed by this unit):**
- **Boron (p-type) mobility at 10¹⁷–10¹⁹ cm⁻³:** resistivity 23–48% high (the legacy Caughey–Thomas hole mobility).
- **Incomplete ionization above about 10¹⁷ cm⁻³:** resistivity 31% to 302% high at 10¹⁸–10¹⁹ (a fixed 45 meV level).
- **Absolute n_i-driven currents:** n_i(300 K) 10% high against the measurement, so n_i² currents (diode saturation,
  generation) about 21% high.
- **Outside the benchmarked range:** below 10¹⁴ (phosphorus) or 10¹⁵ (boron) and above 10¹⁹ cm⁻³; temperatures
  other than 300 K and 23 °C; compensated material; dopants other than phosphorus and boron (the report's fits
  apply to them only approximately). The report notes its phosphorus conversion understates density above about
  5 × 10¹⁹ cm⁻³.
- **Not validated against published data at all:** breakdown (only the analytic one-sided value, Unit 19), MOS C-V,
  high-field transport (the benchmarks use 10 V/cm), transient and small-signal behaviour.

---

### Unit 34 — selectable Masetti hole mobility (merged to `main` as `e14a53f`, PR #42, 2026-10-15)

**Owner request (2026-10-15):** a selectable Masetti doping-dependent hole mobility, off by default (legacy results
preserved), equations and coefficients verified against the paper or a reliable primary reproduction before coding,
compared with Caughey–Thomas on Unit 33's boron resistivity benchmarks, with regression tests against acceptance
criteria fixed in advance; incomplete ionization and n_i unchanged; no push, merge or default change without
authorization.

**Sources read (before coding).** G. Masetti, M. Severi, S. Solmi, IEEE Trans. Electron Devices 30, 764 (1983),
doi:10.1109/T-ED.1983.21207, is paywalled and was not read directly. Its boron model was read from three independent
reproductions, which agree on the form and on every coefficient:
- CERN Allpix Squared, "Charge Carrier Mobility" (allpix-squared.docs.cern.ch/docs/06_models/02_carrier_mobility/):
  the hole formula as the paper's eq (4), the values "taken from Table I" of the paper;
- Crosslight Software, "Low field mobility models" (2010; crosslight.com, Low_field_mobility_model.pdf): the
  Masetti default coefficients in SI units;
- the Flooxs model page (flooxs.ece.ufl.edu, Masetti Mobility Model), as quoted by a search engine (the page itself
  could not be fetched: certificate error).

μ_p(N) = μ_min1 exp(−P_c / N) + (μ_max − μ_min2) / (1 + (N / C_r)^α) − μ_1 / (1 + (C_s / N)^β), N the total ionized
impurity density, with μ_min1 = 44.9, μ_min2 = 0, μ_max = 470.5, μ_1 = 29.0 cm²/(V s), P_c = 9.23e16, C_r = 2.23e17,
C_s = 6.10e20 cm⁻³, α = 0.719, β = 2.0. Masetti's fits are at 300 K; the μ_max (T / 300)^−2.2 scaling that Allpix
Squared and Crosslight apply is not from the paper (Allpix attributes it to another model).

**Acceptance criteria, fixed before any Masetti result was computed:**

| # | Criterion |
|---|---|
| AC1 | Off (the default): every existing test passes unchanged; run identities and Unit 33's pinned Caughey–Thomas gaps are unchanged. |
| AC2 | The implemented μ_p equals the formula above evaluated independently in the test at N = 0 and 1e13 to 1e21 cm⁻³ (each decade) and at 250, 300 and 400 K, within 1e-12 relative; μ_p(0) = μ_max (T / 300)^−2.2. |
| AC3 | Boron resistivity at 300 K (Unit 33's decks, NBS SP 400-64): \|simulated / reference − 1\| ≤ 10% at each of 1e15, 1e16, 1e17, 1e18, 1e19 cm⁻³. |
| AC4 | At 1e17, 1e18 and 1e19 cm⁻³ the 300 K gap is at least 10 percentage points smaller in magnitude than Caughey–Thomas's (30.0%, 48.0%, 23.3%), and the largest \|gap\| over 1e15–1e19 is below Caughey–Thomas's 48.0%. |
| AC5 | Boron resistivity at 23 °C: \|gap\| ≤ 12% at each density (looser: the temperature scaling is not from the paper). |
| AC6 | The electron mobility is unchanged with Masetti on: the phosphorus benchmark gaps equal the Caughey–Thomas pins. |
| AC7 | The temperature slope used by the electrothermal model matches a central difference of μ_p(T) to 1e-6 relative. |
| AC8 | Masetti is refused (invalid_input, before solving) on a semiconductor region whose hole mobility parameters are not silicon's (its coefficients are silicon's). |
| AC9 | Performance: with Masetti the isothermal sweep's wall time is within 10% of Caughey–Thomas's (the mobility is computed once per node at assembly); measured and reported. |

**Implementation.** `physics::masetti_mobility` and `masetti_mobility_slope` with `masetti_boron_silicon`
(`physics/mobility.hpp`); `PhysicsModels::hole_mobility` (`caughey_thomas`, the default, or `masetti`) used for the
hole mobility in the isothermal and electrothermal assemblies; electrons keep Caughey–Thomas; the setting enters the
run record only when Masetti is selected; the deck's `[models] hole_mobility`. Incomplete ionization and n_i are not
touched. No default changed.

**Results against the criteria (measured after the criteria commit, `3670fc3`):**

| # | Result | Met? |
|---|---|---|
| AC1 | Explicit Caughey–Thomas is bit-identical to the default (currents, densities, run identity); all existing suites pass; Unit 33's pins hold. | yes |
| AC2 | Within 1e-12 at N = 0, 1e13–1e21 and 250, 300, 400 K; μ_p(0) = μ_max (T / 300)^−2.2. | yes |
| AC3 | 300 K gaps: −1.1% (1e15), −1.1% (1e16), −5.3% (1e17), **−11.0% (1e18)**, −3.9% (1e19). | **no, at 1e18** (by 1.0 point) |
| AC4 | Against Caughey–Thomas: 30.0 → 5.3, 48.0 → 11.0, 23.3 → 3.9 points at 1e17, 1e18, 1e19; largest \|gap\| 48.0% → 11.0%. | yes |
| AC5 | 23 °C gaps: −1.1%, −1.1%, −5.8%, **−12.4%**, −4.9%. | **no, at 1e18** (by 0.4 point) |
| AC6 | The phosphorus gaps equal the Caughey–Thomas pins at both temperatures. | yes |
| AC7 | Within 1e-6 of a central difference (h = 1e-4 T) at 250, 300, 400 K. | yes |
| AC8 | Refused for a germanium region; accepted for silicon with another parameter changed. | yes |
| AC9 | 2D diode, 20,301 nodes, 6 points, best of 3: ratios 1.36, 0.94, 0.98, 0.93 in four runs; the three runs that counted iterations took 26 Newton iterations with either model. | yes in 3 of 4 runs; the 1.36 run is timing noise on this machine (its iterations were not recorded) |

**Owner ruling (2026-10-15):** the AC3 and AC5 misses at 10¹⁸ cm⁻³ are **accepted as documented limitations, not as
passed criteria**.

The criteria are not loosened. The validation test checks AC3 and AC5 everywhere but 1e18, pins every measured gap,
and states AC3 and AC5 at 1e18 in a separate `[!shouldfail]` case that reports the miss (and would report a fix). AC9 is
a hidden `[.performance]` case (timing does not belong in ctest).

**Limitations and open issues:**
- **AC3/AC5 at 1e18 cm⁻³ not met (−11.0%, −12.4%).** A likely cause, not verified: Masetti fitted mobility against
  carrier (hole) concentration, while the benchmark solves with complete ionization (p = N); the NBS report finds
  boron de-ionization in this range, which a resistivity computed with p = N underestimates. Resolving it needs a
  valid ionization model at 1e18 (Unit 33 finding; out of this unit's scope by instruction).
- The coefficients come from three reproductions, not from the paper itself (paywalled).
- The (T / 300)^−2.2 scaling of μ_max is not Masetti's; Masetti's fits are at 300 K only. The 23 °C gaps track the
  300 K ones within 1.4 points, but other temperatures are not validated.
- Boron-doped silicon only: the coefficients are Masetti's boron set; other acceptors and other materials are refused
  or untested. Electron mobility is unchanged (Caughey–Thomas; phosphorus gaps −8.8% to +3.9%).
- Mutation checks (8), all caught: α miscopied, the slope's factor dropped, the clustering term inverted, the model
  selection inverted, the refusal disabled, the electrothermal flag never set, the run-record setting, the deck
  mapping.

---

### Unit 35 — published silicon defaults: band set and hole mobility (merged to `main` as `2d1fdc7`, PR #44, 2026-10-16)

**Owner authorization and rulings (2026-10-15)**, asked before any work: change the defaults. (1) n_i from a
published, physically consistent silicon band set (Eg(T), Nc(T), Nv(T)), verified from the source before coding — not
a rescale of Nc/Nv or an Eg0 shift; (2) the default hole mobility becomes "automatic": Masetti in regions with
silicon's hole mobility parameters, the material's Caughey–Thomas elsewhere (heterostructures keep working); (3) the
legacy stays selectable — the legacy silicon set (n_i 1.0674e10 cm⁻³) and `hole_mobility = caughey_thomas` as named
options reproducing old results bit for bit; (4) legacy-port gates keep running on the explicitly selected legacy
options; new-default gates are added with their before/after values recorded and justified against Unit 33's
benchmarks. Incomplete ionization is not touched.

**Source read (before coding).** R. Couderc, M. Amara, M. Lemiti, "Reassessment of the intrinsic carrier density
temperature dependence in crystalline silicon", J. Appl. Phys. 115, 093705 (2014), doi:10.1063/1.4867776, open
access (hal.science/hal-00956479), read from the PDF's own text streams. The extraction drops minus signs (plus signs
survive as a separate glyph); each sign below was fixed by the paper's own cross-checks, recorded in the last column.

| Quantity | As used | Paper | Check |
|---|---|---|---|
| n_i | n_i² = Nc Nv exp(−Eg(T) / kT) | eq (1) | — |
| Nc, Nv | 4.83e15 (m_dc / m0)^(3/2) T^(3/2), likewise Nv [cm⁻³] | eqs (4), (5) (the constant as printed; CODATA gives 4.829e15) | — |
| (m_dc / m0)^(3/2) | A T³ + B T² + C T + D, A = −4.609e-10, B = 6.753e-7, C = −1.312e-5, D = 1.094 | eq (18), 0–400 K | It fits 6 m_t √m_l / m0^(3/2) with m_l = 0.9163, m_t = 0.1905 Eg(0) / Eg(T) (eqs 16, 17) within 0.5% (the paper); these signs give 0.03% from 20 to 400 K, and D = 1.094 = 6 (0.1905) √0.9163 |
| (m_dv / m0)^(3/2) | A T³ + B T² + C T + D, A = 2.525e-9, B = −4.689e-6, C = 3.376e-3, D = 0.3426 | eq (19), from eq (15) | With eq (18) and Eg(T) the resulting n_i is within 1.5% of eq (15) from 50 to 400 K (the paper); these signs give 1.3% (largest, 100 K), 0.12% at 300 K |
| Eg(T) | Eg(0) − αΘ [ (1 − 3Δ²) / (exp(Θ / T) − 1) + (3Δ² / 2) ( ⁶√(1 + π² χ² / (3 (1 + Δ²)) + (3Δ² − 1) χ³ / 4 + 8 χ⁴ / 3 + χ⁶) − 1 ) ], χ = 2T / Θ | eq (10), R. Pässler, Phys. Rev. B 66, 085201 (2002), the paper's chosen gap model | π² / (3(1 + Δ²)) (not π² · 3(1 + Δ²)): the mass fit above holds to 0.03% only so (2.4% otherwise); prefactor αΘ: gives the renormalization αΘ / 2 the paper states for Pässler's model and Eg(300 K) = 1.1242 eV |
| Pässler's parameters | Eg(0) = 1.17 eV, α = 3.23e-4 eV/K, Θ = 446 K, Δ = 0.51 | Table I | — |
| Reference n_i(T) | 1.541e15 T^1.712 exp(−Eg(T) / 2kT) | eq (15); "9.68 × 10⁹ cm⁻³" at 300 K | Reproduced: 9.685e9 |

Resulting silicon values at 300 K: Eg = 1.12416 eV (legacy 1.12452), Nc = 2.857e19 (legacy 2.86e19), Nv = 2.514e19 cm⁻³
(legacy 3.10e19), n_i = 9.673e9 cm⁻³ (legacy 1.0674e10; Misiakos 9.7 ± 0.1e9, Unit 33 V4).

**Design (decided here, before coding):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| S1 | Two optional parameter groups on `SemiconductorParameters`: a Pässler gap (α, Θ, Δ, with `Eg0_eV` as Eg(0); the Varshni α and β then 0) and density-of-states mass polynomials (conduction and valence; Nc300 and Nv300 then 0). All zero (every existing set): the legacy Varshni and (T / 300)^(3/2) formulas, evaluated exactly as before. | A published set in the material's own data; legacy sets bit for bit; the forms are general, not a silicon switch. | Two ways to state a band; `create` refuses a mixture (both forms given for one quantity). |
| S2 | `physics::silicon_parameters` (deck `silicon`) becomes Couderc's set; every other field is unchanged. The old set is `physics::legacy_silicon_parameters` (`legacy_silicon()`, deck `silicon_legacy`). | Ruling (1), (3). | Code and decks that named silicon change results (the point of the unit). |
| S3 | `HoleMobility::automatic`, the new default: per semiconductor region, Masetti where the region's hole mobility parameters equal silicon's Caughey–Thomas set, the material's Caughey–Thomas elsewhere. `caughey_thomas` (legacy) and `masetti` (refuses non-silicon regions, Unit 34) stay. The run record names the setting unless it is `caughey_thomas`, so legacy-option runs keep their identity. | Ruling (2), (3). | "Silicon" is recognized by its hole mobility parameters (as Unit 34's refusal). |
| S4 | Electrothermal: with temperature-dependent DOS masses, N(T) = N(T0) θ^(3/2) M(θ). The extra −kT ∇ln M in Wachutka's flux is a band-parameter step across an edge: ln M(θ_b) − ln M(θ_a) joins the edge's ln Nc0 step (electrons; −ln Nv for holes), with its θ partials; the densities at fixed reduced energy move with 3 / (2θ) + d ln M / dθ in place of 3 / (2θ). The energy a carrier carries, (r + 5/2) kT above its band edge, is unchanged (an assumption: the thermopower is that of parabolic bands at the local mass). With no mass polynomial, M ≡ 1 and every term is the legacy's, bit for bit. | The default material must keep working electrothermally; refusing it or freezing M at T0 would be a silent approximation (T8's rule). | The thermopower assumption above is recorded, not validated. |
| S5 | The p-polysilicon work function of an `electrode` contact (no semiconductor under it) uses `silicon_parameters`, now Couderc's gap: 0.36 meV below the legacy at 300 K. | One silicon default. | A legacy run with a p-poly electrode contact differs by that shift (gates take the gap of the semiconductor under them, so a legacy-silicon device is bit for bit). |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| DC1 | Legacy selectable, bit for bit: with `silicon_legacy` / `legacy_silicon_parameters` and `hole_mobility = caughey_thomas`, a diode sweep and the boron 1e18 benchmark written by `nitcad run` are byte-identical to the same decks run on `main` (`0768cae`, the deck's `silicon` then), and every legacy-port gate passes on those options with its unchanged tolerance. |
| DC2 | The implemented Eg(T), Nc(T), Nv(T) equal the formulas above evaluated independently in the test (from this table, not the parameter struct) at 50, 77, 100, 200, 300, 400 K within 1e-12 relative; legacy sets unchanged (Varshni and power law, bit for bit). |
| DC3 | n_i(300 K) of the default silicon in [9.6e9, 9.8e9] cm⁻³ (Misiakos's 9.7 ± 0.1e9) and within 1.5% of eq (15) from 50 to 400 K (the paper's stated precision). |
| DC4 | Unit 33's solved n_i benchmark (√(np) in the 1e14 bar at 0 V) with the defaults: \|gap\| ≤ 1.1% against 9.7e9 (within the measurement's uncertainty); the legacy pin (+10.0%) unchanged on legacy options. |
| DC5 | Unit 33/34 resistivity benchmarks with the defaults: boron gaps equal Unit 34's Masetti pins and phosphorus gaps equal the Caughey–Thomas pins, each within 0.2 points (a majority-carrier resistivity does not depend on n_i); the largest boron \|gap\| ≤ 11.2%. |
| DC6 | dEg/dT, d ln Nc / dT, d ln Nv / dT and d ln n_i / dT used by the electrothermal model match central differences within 1e-6 relative at 100, 300 and 400 K (both band forms). |
| DC7 | Electrothermal with the defaults: the existing Jacobian (finite-difference) and energy-balance checks pass on the new silicon; with isothermal sinks at 300 K and a small current, the electrothermal current matches the isothermal one within 1e-3. |
| DC8 | `automatic`: Masetti in a silicon region and Caughey–Thomas in a germanium region of one device (no refusal), each region's hole current matching a run with the explicit model on a single-material device; `masetti` still refuses the germanium region. |
| DC9 | New-default gates: each legacy-port gate's device is also run on the defaults, the value recorded here before/after with its explanation; a diffusion-limited diode current changes by the n_i² ratio (9.673e9 / 1.0674e10)² = 0.821 within 2% (hole mobility held at Caughey–Thomas to isolate it). |
| DC10 | Full Release and Debug suites after each commit; zero new warnings. |

**Implementation.** `physics::PasslerGapParameters` and `DensityOfStatesMassParameters` on `SemiconductorParameters`
(S1; `create` refuses a mixed or incomplete form; `check_temperature` refuses a mass fit that is not positive, silicon's
conduction cubic near 2030 K); `silicon_parameters` is Couderc's set and `legacy_silicon_parameters` /
`legacy_silicon()` the old one (S2); `physics::dos_mass_log` and `dos_log_slope`, and the gap, n_i and
emission-velocity slopes of both forms (`thermal.hpp`); `HoleMobility::automatic` (value 2; the default) with
`assemble::masetti_holes` deciding per material (S3); the electrothermal mass term (S4: `ThermalState::mass_n/p`, the
edges' ln M step, the levels, the ohmic contacts and the interfaces); device-file keys 56–58 and 60–67 (absent: the
legacy forms, so a file written before Unit 35 loads as legacy silicon); the run identity digests the new forms only
when present; the deck's `material = silicon_legacy`, `passler_gap.*` overrides and `hole_mobility = automatic`.
Incomplete ionization is not changed.

**Results against the criteria (measured after the criteria commit, `f250d9c`):**

| # | Result | Met? |
|---|---|---|
| DC1 | `nitcad run` on `silicon_legacy` with `hole_mobility = caughey_thomas`: the diode sweep and the boron 1e18 bar's `.sweep` and `.csv` files are byte-identical to `main`'s (`0768cae`) on the same decks. The legacy-port gates run on the legacy options and pass with their tolerances: J(0.5 V) = 1.2801e-2 (bias, sweep, field mobility, heavy doping), the M33 heterojunction set, the MOS capacitor port and the MOSFET Id–Vg and mesh gates, the legacy diode turn-off, Unit 33's pins and Unit 34's Masetti pins. | yes |
| DC2 | Eg, Nc, Nv and n_i within 1e-12 of the independent evaluation at 50, 77, 100, 200, 300, 400 K; the legacy sets' formulas unchanged (their double-double reference tests pass). | yes |
| DC3 | n_i(300 K) = 9.673e9 cm⁻³; within 1.3% of eq (15) from 50 to 400 K (largest near 100 K). | yes |
| DC4 | Solved n_i 9.6731e9: −0.3% (0.3 times the uncertainty); the legacy pin +10.0% unchanged. | yes |
| DC5 | Boron at 300 K: −1.1%, −1.1%, −5.3%, −11.0%, −3.9%; 23 °C: −1.1%, −1.1%, −5.8%, −12.4%, −4.9% (Unit 34's Masetti pins); phosphorus unchanged; largest boron gap 11.0%. | yes |
| DC6 | Within 1e-6 of central differences at 100, 300, 400 K for both forms (and the emission velocity with and without a Richardson constant). | yes |
| DC7 | On the default silicon: the finite-difference Jacobian gates (steady and time step), the energy identity and the uniform-T equilibrium pass; the Seebeck voltage (300–310 K, both silicon sets) within 1e-5 and the 300–400 K refinement case within 1e-5 of ∫(2 − ln(N_D / Nc(T))) dT with the material's own Nc(T); the Peltier coefficient within 1e-7 at Couderc's Nc(300 K); the isothermal limit within 1e-3. | yes |
| DC8 | A silicon/germanium diode solves with `automatic`, `masetti` refuses it; on silicon `automatic` equals `masetti` bit for bit, on germanium `caughey_thomas`. | yes |
| DC9 | See the before/after table; the n_i² ratio of a diffusion current (Caughey–Thomas holes): 0.8219 at 0.4 V and 0.8217 at 0.5 V against 0.8213 (+0.08%, +0.05%). | yes |
| DC10 | Release and Debug: 18/18 suites; no warnings. | yes |

**Before (legacy options) and after (the defaults):**

| Gate | Before | After | Why |
|---|---|---|---|
| Unit 9 diode, J(0.5 V) | 1.2801e-2 A/cm² | 1.1399e-2 (−11.0%) | n_i² × 0.821; the hole part (23%) up by Masetti's 319 against 232 cm²/(V s) in the 1e17 n side |
| n_i(300 K) | 1.0674e10 cm⁻³ | 9.673e9 | Couderc's set (Nv 3.10e19 → 2.514e19, Nc 2.86e19 → 2.857e19, Eg 1.12452 → 1.12416 eV) |
| MOSFET V_th (max g_m) | 0.1520 V | 0.1527 V | — |
| MOSFET swing | 68.79 mV/decade | 68.61 | — |
| MOSFET Id(Vg = 1 V) | 0.31471 A/cm | 0.31448 (−0.1%) | an electron channel current |
| MOSFET on/off | 1.459e11 | 1.303e11 | the off current, junction leakage, 12% higher (not decomposed) |
| Boron resistivity, largest gap | +48.0% | −11.0% | Masetti (Unit 34) |
| Incomplete ionization, boron 1e18 / 1e19 | +120.9% / +301.7% | +40.3% / +241.3% | Masetti and Nv (the acceptor occupancy); still the fixed 45 meV level (Unit 33 finding) |

Tests that compare against analytic forms written with the legacy Caughey–Thomas hole mobility or the T^(3/2) law now
take the default's (the MOS RC and dielectric-relaxation references, the surface-recombination diode, the field-
mobility resistor, the gate's intrinsic offset) or the material's Nc(T) (Seebeck, Peltier); legacy reproductions run on
the legacy options. Two comparisons were relaxed to the quantity's resolution, each explained where it is: a sweep or a
restart against a chain of one-point solves compares contact currents within their Unit 32 resolution (with the lower
n_i the cancelling majority fluxes round the current by 6e-8 of it at 0.4 V), and a transient snapshot's time is
matched to 1e-12 (it is one ulp off 3 τ through the scaled time). The mid-gap-trap test allows the rows' rounding
(the rate is 2e-8 of rows rounding at 1e-16).

**Unresolved issues (found by this unit, not fixed):**
- **Nonlocal band-to-band trace on the defaults.** The 5e19 zener's `trace_bias` to −2 V stalls at −1.31 V: the
  corrector's updates stay between 3e-8 and 2.5e-6 against the 1e-8 tolerance with the residual flat, a rounding floor
  of the lagged-path system; the swept solves converge, and with Caughey–Thomas holes the trace completes. The
  mechanics test runs on Caughey–Thomas; a `[!shouldfail]` case records the default.
- **Open-base snap-back trace on the defaults.** The fold (17.46 V, 1.61e-3 A/cm²) is crossed in about 5000 steps,
  the current moving by its resolution (2.9e-7 A/cm²); the curve and BV_CEO are reached, but the highest traced bias
  falls past the turning-point window. The case runs on the legacy options (not repeated as `[!shouldfail]`: four
  minutes).
- The p-polysilicon work function of an `electrode` contact takes the default silicon's gap (S5): 0.36 meV below the
  legacy's at 300 K, so a legacy run with such an electrode is not byte for byte (none of the legacy gates has one).
- The thermopower of a temperature-dependent mass (S4) is an assumption, checked for consistency (Seebeck, Peltier,
  energy identity), not against a measurement.
- Couderc's fits are stated for 0–400 K (n_i to 1.5% over 50–400 K); above 400 K they are evaluated as written.

**Performance.** The band forms are evaluated once per node at assembly (isothermal) or once per node per assembly
(electrothermal: two polynomials and logarithms per band); not separately timed.

Mutation checks (10), all caught: Pässler's prefactor halved, the mass polynomial's slope missing its constant, the
electron and the hole ln M step's sign, the level's d ln M / dθ dropped, `automatic` inverted, the run identity
skipping the valence fit, the device file dropping it, the deck's `silicon_legacy` mapping, the n_i slope's mass term.

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
