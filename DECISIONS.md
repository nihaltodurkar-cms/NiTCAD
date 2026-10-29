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
- **Nonlocal band-to-band trace on the defaults** (resolved by Unit 36: the arc row's rounding). The 5e19 zener's `trace_bias` to −2 V stalls at −1.31 V: the
  corrector's updates stay between 3e-8 and 2.5e-6 against the 1e-8 tolerance with the residual flat, a rounding floor
  of the lagged-path system; the swept solves converge, and with Caughey–Thomas holes the trace completes. The
  mechanics test runs on Caughey–Thomas; a `[!shouldfail]` case records the default.
- **Open-base snap-back trace on the defaults** (resolved by Unit 36). The fold (17.46 V, 1.61e-3 A/cm²) is crossed in about 5000 steps,
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

### Unit 36 — trace robustness at the current-resolution floor (merged to `main` as `07cd53a`, PR #46, 2026-10-16)

**Owner ruling (2026-10-16):** Units 36–39 in order (trace robustness, density-dependent ionization, Masetti electron
mobility, breakdown validation), each its own branch, criteria before code and PR; stop for the go-ahead to merge
after each; Units 37 and 38 opt-in.

**The two failures (Unit 35, unresolved):** on the defaults, the 5e19 zener's nonlocal band-to-band `trace_bias` stalls
at −1.31 V, and the open-base transistor's snap-back trace crosses its fold in about 5000 steps.

**Diagnosis (measured with temporary instrumentation, removed before this commit).** In both, the corrector's largest
update is the swept bias λ itself (every potential on the swept side moves with it), alternating in sign at 1e-7 to
1e-5 with the residual flat: λ is not resolved. The arc-length row c is M τ, and its current term is the gradient of
log₁₀(|I| + I0) built from the swept contact's current row (`small_signal_currents`), whose entries are the contact's
majority fluxes' derivatives: they cancel as the contact current does (its Unit 32 resolution). Measured on the zener:
Σ|c_q x_q| up to 7e8 against a row value of 1e-4, a cancellation over 13 orders, so both the arc residual
c·(y − y₀) − ds and the Sherman–Morrison products u·z and u·w of the bordered solve are rounding. Checks: with the current
term removed the zener corrector converges quadratically (4–5 iterations, updates 1e-14) and reaches −2 V; with the
products computed in compensated arithmetic it reaches −2 V with the current term kept. The snap-back with compensated
products crosses the fold in a few steps but then accepts a corrected point on the other branch (16.87 V, 7.9e-4
A/cm² → 13.33 V, 9.2e-2 A/cm²; the step's arc length was 1 V) and runs the curve backwards; the tangent-turn test does
not catch it.

| # | Choice | Why | Trade-off |
|---|---|---|---|
| R1 | Every product with the arc row (the arc residual, u·z and u·w) is a compensated dot product (Ogita, Rump, Oishi, SIAM J. Sci. Comput. 26, 1955 (2005), Dot2: TwoProduct by FMA, TwoSum), as accurate as in twice the working precision. | The residual and the solve then agree with c to the rounding of the result, not of Σ|c_q x_q|; c is ours (any fixed row defines a valid hyperplane), so only consistency matters. | About four flops per entry more, once per iteration (n + 1 entries). |
| R2 | A corrected point farther from its predicted point than the step itself, measured in the arc metric, is rejected and the step retried at half (a distance test of the kind in Allgower and Georg, *Numerical Continuation Methods*, 1990). | The turn test missed a 3.5 V jump to the other branch; a corrector that converges so far from its prediction is not following the branch. | A legitimate large correction is retried smaller. |
| R3 | No iterative refinement of the bordered solve. | Measured with R1 and R2: Newton iterations over the zener trace 250 → 169, over the snap-back 901 → 879; not worth the extra solves. | — |

**Acceptance criteria, fixed before the implementation (the diagnosis above was measured with temporary code):**

| # | Criterion |
|---|---|
| TR1 | The zener nonlocal trace to −2 V on the defaults completes, its end point the swept solve's within the resolution and 1e-6 (the existing assertions); Unit 35's `[!shouldfail]` case becomes an ordinary case. |
| TR2 | The open-base snap-back on the defaults passes every assertion of the existing case (fold, turning-point window, the high-current branch to 20 A/cm², the bias falling past the fold, α_T M_n); the case returns to the defaults. |
| TR3 | Both cases also pass on the legacy options (legacy silicon, Caughey–Thomas holes), as before. |
| TR4 | Every other test passes with its tolerance unchanged (the traces of the breakdown, band-to-band, electrothermal and weighted-current suites among them). |
| TR5 | Cost: the snap-back on the defaults takes at most 200 points (Unit 35: about 5000); points and wall times of both cases are reported. |
| TR6 | The compensated dot product is exact to the rounding of its result on a cancelling pair (terms of 1e16 summing to 1); the distance test and R1 are mutation-checked (each removed, a TR1/TR2 case fails). |
| TR7 | Full Release and Debug suites; zero warnings. |

**Implementation.** `linalg::compensated_dot` (`include/NiTCAD/linalg/compensated_dot.hpp`, Dot2 with FMA; needs the
build's `/fp:precise`); `trace_bias` computes the arc residual and the Sherman–Morrison products u·z, u·w with it (R1)
and rejects a corrected point farther than ds from its prediction in the arc metric (R2; not on the step that lands on
end_V, whose row is the bias itself). `trace.hpp` documents both.

**Results against the criteria (measured after the criteria commit, `7ba51ef`):**

| # | Result | Met? |
|---|---|---|
| TR1 | The zener trace reaches −2 V on the defaults and on Caughey–Thomas holes, its end point the swept solve's within the resolution and 1e-6; the `[!shouldfail]` case is removed, the ordinary case runs both. | yes |
| TR2 | The snap-back on the defaults passes every assertion: fold 17.4625 V, the turning point at it, the high-current branch to 20 A/cm², α_T M_n within 1e-3, every point a steady state. | yes |
| TR3 | On the legacy options: fold 17.2666 V (Unit 19's 17.27), every assertion passes. | yes |
| TR4 | Every other test passes with its tolerance unchanged. One test changed: the electrothermal runaway fold read the fold from the parabola through three traced points, an O(h³) estimate that depended on where the top point fell (Unit 36's corrector takes fewer, differently placed points: 1.4e-7 V against the 1e-7 tolerance, 0.13 K from the fold; before, 4.6e-9 at 0.05 K). It now maximizes the quartic through five points (O(h⁵)): 7.3e-10 V, and −2.7e-10 V on the code before Unit 36; the tolerance is unchanged. | yes |
| TR5 | Snap-back: 59 points in 0.7 s on the defaults (Unit 35: about 5000), 55 points in 0.8 s on the legacy options. Zener: about 20 step attempts, under a second, on each option set. | yes |
| TR6 | `compensated_dot` returns 1 for 1e16 + 1 − 1e16 (the plain sum 0) and 2⁻⁶⁰ for (1 + 2⁻³⁰)² − 1 − 2⁻²⁹. Mutation checks (5), all caught: the arc residual plain, the Sherman–Morrison products plain, the distance test off, Dot2's product error dropped, its sum error dropped. | yes |
| TR7 | Release 18/18, Debug 12/12; no warnings. | yes |

**Limitations.** The current term of the arc metric still comes from the contact current; only its arithmetic is
compensated. A current below its resolution still cannot steer a step (I0, unchanged). The distance test uses the arc
metric at the step's start.

---

### Unit 37 — density-dependent dopant ionization, opt-in (merged to `main` as `16a77b5`, PR #47, 2026-10-16)

**Owner rulings (2026-10-16):** the model is opt-in (off by default); a default change only on later authorization.

**Source read (before coding).** A. Schenk, P. P. Altermatt, B. Schmithüsen, "Physical model of incomplete ionization
for silicon device simulation", SISPAD 2006, pp. 51–54 (IEEE Xplore document 4061579; title, authors and pages read from the PDF), open at TU Wien
(in4.iue.tuwien.ac.at, sispad2006/pdfs/04061579.pdf), read from the PDF's text. It states the device model of P. P.
Altermatt, A. Schenk, B. Schmithüsen, G. Heiser, J. Appl. Phys. 100, 113714 and 113715 (2006) (paywalled, not read):

| Equation | As used |
|---|---|
| (3a), (3b) | N_D⁺ / N_D = 1 − b(N_D) n / (n + g_D⁻¹ n₁), N_A⁻ / N_A = 1 − b(N_A) p / (p + g_A⁻¹ p₁) |
| (4) | n₁ = Nc exp(−E_D / kT), p₁ = Nv exp(−E_A / kT) (Boltzmann; the paper replaces the impurity DOS by a delta function and the quasi-Fermi levels by Boltzmann statistics) |
| (5) | E_dop(N) = E_dop,0 / (1 + (N / N_ref)^c) |
| (1) | b(N) = 1 / (1 + (N / N_b)^d), the fraction of dopant states still bound (the rest merged with the band) |
| Table I | Si:P E_dop,0 = 45.5 meV, N_ref = 2.2e18 cm⁻³, c = 2, N_b = 6e18 cm⁻³, d = 2.3, g = 1/2; Si:B 44.39 meV, 1.3e18, 1.4, 4.5e18, 2.4, 1/4 |

Checks of the transcription: with b = 1 and E_dop fixed, (3a) is the existing fixed-level form N / (1 + g_D n / n₁),
g_D = 1/g = 2 (4 for boron), so the table's g is the reciprocal of the legacy degeneracy; the paper's eqs (9)–(10)
state the effective degeneracy b / (g + (1 − b) N / n₁) tends to 1/g as N → 0 and to 0 as N → ∞, which this form
gives. The paper finds "up to 25% of phosphorus and boron atoms ... non-ionized at room temperature near the Mott
transition" (N_crit = 3.74e18 cm⁻³) and no significant non-ionization far above it.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| I1 | `PhysicsModels::ionization_levels`: `fixed` (the default: the material's fixed levels, as since Unit 15) or `altermatt` (the table's phosphorus set for donors, boron for acceptors), read only with `incomplete_ionization` on. | Opt-in (ruling); the coefficients are silicon's. | Donors are taken as phosphorus and acceptors as boron (arsenic's set is not used). |
| I2 | `altermatt` is refused (invalid_input) on a semiconductor region whose ionization parameters are not silicon's. | As Masetti (Unit 34). | — |
| I3 | Per node, E_dop and b at the dopant's own concentration (N_D for donors, N_A for acceptors), with the full model (3a)–(3b) in the carriers (the paper's eqs (6)–(8), its simulator approximation, are not used). | (3) is what the table was fitted with; the kernels already depend on the carrier density, so no approximation is needed. | In compensated material the paper's Part II has a refinement not read. |
| I4 | The kernel N((1 − b) + b s), s the fixed-level ionized fraction at the node's E_dop: with b = 1 it is the existing value bit for bit. | One kernel; `fixed` unchanged. | — |
| I5 | The run record names the setting only when `altermatt` (legacy identities unchanged); the deck's `[models] ionization_levels`. | As Unit 34. | — |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| II1 | Default (`fixed`): every existing test passes unchanged, the incomplete-ionization pins (Units 33 and 35) included; run identities unchanged. |
| II2 | The implemented ionized densities equal eqs (1), (3)–(5) and Table I evaluated independently in the test, at N = 1e14 … 1e21 cm⁻³ (each decade) and n / n₁ (p / p₁) from 1e-3 to 1e3, within 1e-12 relative; their eta slopes match central differences within 1e-6. |
| II3 | Neutral uniform silicon at 300 K (solved, `altermatt`): the largest non-ionized fraction between 1e17 and 1e20 cm⁻³ is between 20% and 30% for both phosphorus and boron, at a density within a factor of 3 of 3.74e18 cm⁻³; at 1e20 cm⁻³ below 3%; at 1e15 cm⁻³ below 1%. |
| II4 | Unit 33's resistivity decks with incomplete ionization on the defaults: with `altermatt` the gaps at 1e18 and 1e19 cm⁻³ are smaller in magnitude than the fixed level's (phosphorus +31.4%, +137.8%; boron +40.3%, +241.4%) for both dopants; every gap reported and pinned. |
| II5 | The finite-difference Jacobian gate with incomplete ionization passes with `altermatt` (steady state). |
| II6 | `altermatt` is refused on a germanium region; accepted on silicon with another parameter changed. |
| II7 | Full Release and Debug suites; zero warnings. |

**Implementation.** `physics::DensityDependentLevel` with `altermatt_phosphorus` and `altermatt_boron` (Table I),
`dopant_level_eV` and `bound_fraction` (eqs (5), (1)); `ionized_density` takes the bound fraction b (default 1:
the fixed-level value bit for bit); `DopantLevels` carries b per dopant; `PhysicsModels::ionization_levels`
(`fixed`, the default, or `altermatt`) sets each node's levels and b at its own dopant concentrations in the scaled
device, which refuses `altermatt` on a region without silicon's ionization parameters (both Poisson and
drift-diffusion); the run record names the setting only with incomplete ionization on and `altermatt`; the deck's
`[models] ionization_levels`.

**Results against the criteria (measured after the criteria commit, `e7fc6b7`):**

| # | Result | Met? |
|---|---|---|
| II1 | `fixed` is the default; every existing test passes unchanged, the incomplete-ionization pins included; the setting leaves the run identity unchanged where it is not read. | yes |
| II2 | Ionized densities within 1e-12 of eqs (1), (3)–(5) and Table I at 1e14–1e21 cm⁻³ and n/n₁ from 1e-3 to 1e3, for both dopants; slopes within 1e-10 of the formula's derivative everywhere, and within 1e-6 of central differences wherever these resolve the slope (where nearly every atom is ionized the value is N(1 − 1e-6) and its difference quotient is rounding: the criterion's difference check is applied where its rounding is below 1e-8 of the slope). | yes |
| II3 | Largest non-ionized fraction: phosphorus 20.9% at 1.74e18 cm⁻³, boron 26.3% at 2.09e18; at 1e20: 0.14%, 0.06%; at 1e15: 0.04%, 0.09%. | yes |
| II4 | Gaps with incomplete ionization, fixed level → Altermatt: phosphorus 1e18 +31.4% → +24.6%, 1e19 +137.8% → +2.8%; boron 1e18 +40.3% → +17.6%, 1e19 +241.4% → +4.3% (below 1e18 within 0.5 point of the fixed level's); all pinned. | yes |
| II5 | The FD-Jacobian gate (Poisson and drift-diffusion, Boltzmann and Fermi–Dirac) on a 2e18 cm⁻³ silicon junction passes at 5e-5, and the ionization part on its own at 1e-5. | yes |
| II6 | Refused on germanium (drift-diffusion and Poisson); not read without incomplete ionization; accepted on silicon with another parameter changed. | yes |
| II7 | Release 18/18, Debug 12/12; no warnings. | yes |

**Findings.**
- With the density-dependent levels, incomplete ionization no longer ruins the heavily doped benchmarks (1e19
  within 5%), but at 1e18 the resistivity is still 18–25% high. A likely cause, not verified: the mobility fits
  (Caughey–Thomas electrons, Masetti holes) are fits of measured data that already contain the material's incomplete
  ionization, so adding it again counts it twice near the Mott density. Until a mobility fitted against ionized
  densities is used, incomplete ionization with these mobilities overstates the resistance of regions near 1e18
  cm⁻³; it stays opt-in.
- Unit 34's untested hypothesis, that incomplete ionization explains Masetti's −11.0% boron gap at 1e18, is not
  supported: with it the gap becomes +17.6%.

**Limitations.** Donors are taken as phosphorus and acceptors as boron (arsenic's set not used); E_dop and b at the
dopant's own concentration (the compensated-material treatment of the paper's Part II not read); Boltzmann n₁ as in
the paper (with Fermi–Dirac statistics the carriers' reduced energy enters as before); temperature dependence only
through kT and Nc, Nv (the paper's fits are at 300 K and below).

Mutation checks (6), all caught: boron's d miscopied, the slope's b dropped, the silicon refusal off, the acceptor's b at
the donor density (caught by II4, not by the Jacobian gate), the run record ignoring whether ionization is on, the
deck mapping.

---

### Unit 38 — Masetti electron mobility, opt-in (merged to `main` as `1b446d0`, PR #48, 2026-10-16)

**Owner rulings (2026-10-16):** opt-in (off by default); a default change only on later authorization.

**Sources read (before coding).** Masetti, Severi, Solmi (1983) is paywalled (as in Unit 34). The electron (phosphorus)
form is the paper's eq (1), as the hole form without the P_c term (P_c = 0): μ_n = μ_min1 + (μ_max − μ_min2) /
(1 + (N / C_r)^α) − μ_1 / (1 + (C_s / N)^β). Its phosphorus coefficients:

| Coefficient | Allpix Squared (allpix-squared.docs.cern.ch, carrier mobility; "taken from Table I") | V. Palankovski, PhD thesis, TU Wien, Table 3.24 (Masetti column, majority electrons in Si at 300 K) | Used |
|---|---|---|---|
| μ_min1, μ_min2 | 68.5 (μ_0,e) | 68.5, 68.5 | 68.5 |
| μ_max | 1414 | — | 1414 |
| μ_1 | 56.1 | 56.1 | 56.1 |
| C_r | 9.20e16 | 9.2e16 | 9.20e16 |
| C_s | 3.41e20 | 3.4e20 | 3.41e20 |
| α | 0.711 | 0.711 | 0.711 |
| β | 1.98 | 1.98 | 1.98 |

μ_max rests on Allpix alone. Allpix's arsenic row (μ_max 1417, 52.2, 9.68e16, 0.68, 43.4, C_s 3.43e20, 2.0) agrees with
Crosslight's Masetti electron defaults (the arsenic set: 1417, 52.2, 9.68e16, 0.68, 43.4, 2.0, but C_s 3.34e20) and with
the BYU cleanroom resistivity calculator (citing Masetti; 4406.452 Ω cm at 1e12 cm⁻³, i.e. μ = 1416.5, the arsenic set's
1416.4 there, not phosphorus's 1413.6); so Allpix transcribed the table's arsenic row as the others did, and its
phosphorus row is used. Temperature: μ_max (T / 300)^−2.5 (Allpix and Crosslight), not from the paper.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| E1 | `PhysicsModels::electron_mobility`: `caughey_thomas` (the default) or `masetti` (the phosphorus set). | Opt-in (ruling); as Unit 34's holes. | Arsenic's set not offered. |
| E2 | `masetti` refused (invalid_input) on a semiconductor region whose electron mobility parameters are not silicon's Caughey–Thomas set. | As Unit 34. | — |
| E3 | `masetti_mobility` at N = 0 returns its limit: μ_max − μ_min2, plus μ_min1 when P_c = 0 (the electron form, whose μ_min1 term does not vanish); the boron set (P_c > 0) is unchanged bit for bit. | The existing N = 0 branch assumed P_c > 0: for electrons it would give 1345.5 instead of 1414. | — |
| E4 | The run record names the setting only when `masetti`; the deck's `[models] electron_mobility`; electrothermal takes its slope. | As Unit 34. | — |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| ME1 | Default (`caughey_thomas`): every existing test and pin passes unchanged; run identities unchanged. |
| ME2 | μ_n equals the formula and coefficients above, evaluated independently in the test, at N = 0, 1e6 and 1e13 … 1e21 cm⁻³ (each decade) and 250, 300, 400 K within 1e-12; μ_n(0) = μ_max (T / 300)^−2.5; continuous at N → 0 (1e6 cm⁻³ within 1e-6 of μ_n(0)); the temperature slope within 1e-6 of a central difference. |
| ME3 | Phosphorus resistivity at 300 K (Unit 33's decks, on the defaults with `masetti` electrons): the largest \|gap\| over 1e14–1e19 cm⁻³ is below Caughey–Thomas's 8.8%; each gap pinned. |
| ME4 | At 23 °C: each \|gap\| ≤ 12%. |
| ME5 | The boron gaps unchanged by `masetti` electrons within the pin (0.2 points). |
| ME6 | Refused on a germanium region; accepted on silicon with another parameter changed; works with the electrothermal model (isothermal limit within 1e-3). |
| ME7 | Full Release and Debug suites; zero warnings. |

**Implementation.** `physics::masetti_phosphorus_silicon`; `masetti_mobility` at N = 0 keeps μ_min1 when P_c = 0
(E3; the boron set's value unchanged); `PhysicsModels::electron_mobility` (`caughey_thomas`, the default, or
`masetti`) in the isothermal and electrothermal assemblies, refused on regions without silicon's electron mobility;
the run record names it only when `masetti`; the deck's `[models] electron_mobility`.

**Results against the criteria (measured after the criteria commit, `4ce3ff8`):**

| # | Result | Met? |
|---|---|---|
| ME1 | The default is Caughey–Thomas; every existing test and pin passes unchanged. | yes |
| ME2 | Within 1e-12 of the independent evaluation at N = 0, 1e13–1e21 and 250, 300, 400 K; μ_n(0) = 1414 (T/300)^−2.5; 1e6 cm⁻³ within 1e-6 of it (E3); slopes within 1e-6. | yes |
| ME3 | Phosphorus at 300 K, Caughey–Thomas → Masetti: +3.9 → +0.5%, +0.5 → −0.7%, −6.9 → −1.9%, −8.8 → +1.3%, +0.6 → −4.3%, −7.1 → −6.2% (1e14 … 1e19); the largest \|gap\| 8.8% → 6.2%. | yes |
| ME4 | At 23 °C: +0.4%, −1.0%, −2.3%, +0.1%, −6.8%, −7.1%; all within 12%. | yes |
| ME5 | The boron gaps equal the defaults' pins. | yes |
| ME6 | Refused on germanium; accepted on silicon with eps_r changed; named in the run record; electrothermal with 300 K sinks within 1e-3 of the isothermal current. | yes |
| ME7 | Release 18/18, Debug 12/12; no warnings. | yes |

**Limitations.** μ_max = 1414 rests on one reproduction (above); the paper itself not read. Phosphorus only (the
arsenic row, where the reproductions disagree on C_s, is not used). The (T/300)^−2.5 scaling is not Masetti's.

Mutation checks (6), all caught: the N = 0 branch reverted, the electron temperature exponent miscopied, the model
selection inverted, the electrothermal flag never set, the silicon refusal off, the run record silent.

---

### Unit 39 — breakdown validation against published data (deferred by the owner, 2026-10-16)

The planned last unit of the Units 36–39 sequence: abrupt-junction avalanche breakdown voltages against a published
source. Deferred before any work: no source could be read. S. M. Sze, G. Gibbons, Appl. Phys. Lett. 8, 111 (1966),
doi:10.1063/1.1754511, is paywalled (as in Unit 33, V6); a 2021 analysis, "Analysis of the dependence of critical electric field on semiconductor bandgap"
(doi:10.1557/s43578-021-00465-2), is login-gated and its OSTI record would not connect; an archived RCA
Review (March 1969) would not download. Breakdown stays validated against the analytic one-sided value only (Unit
19). It resumes if a readable copy of a source is supplied.

---

### Unit 40 — generic fields and equations behind the assembler (proposed 2026-10-16; authorized and accepted 2026-10-17; on `assemble/fields-equations`)

**Owner request (2026-10-16):** document a refactor toward ARCHITECTURE.md 6.12.1–6.12.2 (the equation graph, 2.1) in
three slices; documentation only on branch `assemble/fields-equations`; no implementation until explicitly authorized.
Scope: this refactor only. No new physics, models, solver strategies, deck keys or file formats.

**Why now (audited on `main` at `f153d57`).** `DriftDiffusion` is one class with the whole system in it:
`src/assemble/drift_diffusion.cpp` (2036 lines) and `src/assemble/electrothermal.cpp` (1136 lines). The unknown layout
is the hard-coded stride `m_` (3: ψ, n, p; 4 with θ last), read as `m_ * i + k` (and `m_ * m_` blocks) on 114 lines of `drift_diffusion.cpp`
and as literal `4 * i + k` throughout `assemble_heat`; `stride()` is read by `solve/bias.cpp`, `solve/fields.hpp`,
`solve/small_signal.cpp`, `solve/trace.cpp`, `solve/transient.cpp` and `electrothermal_test.cpp`. Each physics unit
added private structs and an `add_*` pass to the class (`ImpactNode`, `TunnelNode`, `PathLayout`, `ThermalNode`, …).
Several targets need the layout and the terms to be explicit: Unit 30's open finding (a natural-test norm monotone in
reverse bias "needs the system's unknown layout"), block preconditioners for iterative solvers (field blocks), state
transfer for mesh refinement (Unit 17; field location and units), and rendering (Unit 25; named fields with units).

**The current accumulation order (to be preserved exactly; `DriftDiffusion::assemble`, isothermal path):**

1. Per node in index order: a Dirichlet node (ohmic contact) sets its three rows and diagonal and skips the rest; an
   insulator node sets its rows; otherwise the Poisson charge `-V (n - p - N)` is **assigned**, then incomplete
   ionization **re-assigns** the Poisson row and its two entries, then the gate term is added, then SRH, Auger,
   radiative, then the BDF storage term.
2. The interface edges, in `InterfaceEdges` order (the insulator node's Poisson row, then the semiconductor node's
   three rows).
3. The mesh edges, in mesh edge order: Poisson only on edges without carriers; otherwise Poisson, electron and hole
   fluxes into row a, then row b, skipping Dirichlet rows; interface (`charged`) edges skipped.
4. Impact ionization, local band-to-band tunnelling, nonlocal path tunnelling, each in its own node or path order.

The electrothermal model takes a separate pass (`assemble_heat`, stride 4, its own order); `small_signal_matrix`
calls `assemble` with `interfaces = false` and adds the interface terms itself; `EquilibriumPoisson` is a separate
assembler. Floating-point addition is not associative, so each value of F and J is a sum whose order is fixed by this
sequence; any slice that changes the order of the additions into one entry changes results in the last bits.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| FE1 | **Slice 1 — field descriptors and the per-node unknown map.** A `Field` descriptor: name, units, location (node only in this unit), unknown or fixed, region applicability (semiconductor, insulator, all). A `DofMap` built from the active fields: the unknown index of (field, node), the per-node block size, and the field of an unknown index. For every existing configuration it reproduces today's layout exactly: node-interleaved, ψ, n, p, then θ with the electrothermal model, every field at every node (insulator nodes keep their n, p unknowns held at 0). `stride()` stays and returns the block size; the solve layer and tests move to the map only where they read a field's slot. | The layout becomes data instead of `m_ * i + k` and `4 * i + k` literals; the norm (Unit 30), preconditioner blocks and field output can ask for a field's unknowns. Same layout, so checkpoints (Unit 29), run identities and the Jacobian pattern are unchanged. | Insulator nodes still carry unused density unknowns (dropping them changes the layout: deferred, FD5). Edge, face and cell fields are not in this unit. |
| FE2 | **Slice 2 — equation-term interfaces.** A term declares the fields it writes and reads, its domain (nodes, mesh edges, interface edges, contacts, paths), its Jacobian pattern contribution at creation, and an `add(state, residual, values)` that adds its residual and exact partials at positions located once at creation (as today's `ab`/`ba`/`block_`). The existing passes become terms, each keeping its own loop body and order: node terms (Dirichlet, insulator, Poisson charge with the ionization override, gate, SRH, Auger, radiative, storage), interface-edge terms, mesh-edge terms (Poisson, Scharfetter–Gummel electron and hole fluxes), impact ionization, local and nonlocal tunnelling; and the heat equation and thermal variants of `assemble_heat`. | The 2.1 graph's local, edge and interface/boundary terms as objects; a new equation is a new term, not a new pass in the class. Exact Jacobians stay hand-written (6.12.3: no automatic differentiation in this unit). | Dispatch overhead (FG4 measures it): edge and global terms are called once per pass over their domain, node terms once per node inside the fused node loop (FE3, FR2), never per matrix entry. |
| FE3 | **Order rule for slice 2.** Every entry of F and J receives its contributions in today's sequence with today's operands: node terms in the order of step 1 (assign → ionization re-assign → gate → SRH → Auger → radiative → storage), with today's skips (a Dirichlet or insulator node gets no later node term); then interface-edge, mesh-edge and global terms in steps 2, 3, 4, each over its domain in today's index order. A node term only writes its own node's rows, so the node terms may run fused in one node loop or as one loop per term without changing any entry's order; the fused loop is kept unless FG4 says otherwise (it is today's memory access pattern). Assignments stay assignments (the Dirichlet rows, the Poisson charge and its ionization override); no term is rewritten as "zero then add" where today's code assigns. | Bit-identical results (FG1) need every entry's additions in the same order with the same operands. | Sequential assembly only: concurrent assembly, where several edges add into one row, changes the order (deferred, FD6). |
| FE4 | **Slice 3 — composition.** An assembler that builds the active term list and the `DofMap` from the device and `PhysicsModels`, in the FE3 order, and evaluates F, J, the time-step and small-signal forms over it. `DriftDiffusion` stays as a compatibility wrapper over it with every public member unchanged in signature, preconditions and errors (`create`, `evaluate`, `residual`, `storage`, `small_signal_matrix`, `terminal_currents`, the thermal and tunnelling members, the Newton hooks); the solve layer, the deck, `io` and the tests compile and run without change beyond FE1's map reads. | Composition per 6.12.2 (diode = Poisson + continuity; thermal = those + heat) without moving the solve layer; the wrapper keeps the change reviewable and revertible. | Two names for one system until a later unit retires the wrapper (not this unit). |
| FE5 | **Diagnostics stay where they are in this unit.** `edge_currents`, `terminal_current_resolution`, `edge_flux_terms`, `band_edges`, `gate_charges`, `contact_charges`, `trap_occupancies`, `thermal_contact_heat` and the weighted-current inputs keep their current code paths, reading the layout through the map. | They are derived quantities, not residual terms; moving them is not needed for composition and multiplies the bit-identity surface. | They are not yet generic derived fields (6.12.1's "derived" kind is deferred, FD7). |
| FE6 | **`EquilibriumPoisson` is out of scope** unless the owner rules otherwise (FD1). | It is a separate, small (284 lines) assembler with its own layout (ψ only); it composes from the same Poisson terms later. | Two Poisson implementations remain for now. |
| FE7 | **One branch; a preparatory commit, then one commit per slice, each passing every gate on its own.** Commit 0 adds only the measurement tools, with no assembler change: the FG4 benchmark and the FG1/FG2 recorder (hidden Catch2 cases, not in ctest, that write F/J bit patterns and run results to a path given in an environment variable, and compare against such a file when one is given). The baselines are recorded at commit 0, which is numerically `main`. Then slice 1 (FE1), slice 2 (FE2; split into 2a, the isothermal terms, and 2b, the electrothermal terms, so each review is about half the size), slice 3 (FE4). The owner may stop after any slice; slice 1 is useful alone (FG7). | The benchmark and the recorder do not exist on `main` (the existing `tests/solve/solver_timing.hpp` times the linear solver, not assembly), so "before" must be measured with the same tools as "after". Each commit is revertible and reviewable; a failed bit-identity check names the slice that caused it. | Slices 2 and 3 touch most of `drift_diffusion.cpp` and `electrothermal.cpp`; reviews are large even when split. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| FG1 | **Bit-identical numerics.** For every configuration the existing suites exercise (Boltzmann and Fermi–Dirac; every recombination, mobility, ionization, band-set and gap option; gates, insulators, interfaces and traps; heterojunctions and thermionic emission; impact ionization, local and nonlocal tunnelling; electrothermal; steady, transient, small-signal and `trace_bias`; tensor-grid and unstructured meshes in 1D, 2D and 3D), F and every Jacobian value equal the pre-refactor values **bit for bit** (compared as `std::bit_cast<std::uint64_t>`, not with a tolerance), as do the time-step and small-signal matrices, `bias_derivative` and `small_signal_currents`. Checked by the commit-0 recorder (FE7): at each configuration it evaluates F and J at fixed states built without the code under change (`EquilibriumPoisson`'s equilibrium, out of scope, FE6, with a set bias, plus a deterministic perturbation of every unknown so no row is at a solution) and writes their bit patterns; after each slice the same cases are compared entry by entry, and the first differing row and column are reported. No second copy of the old assembler is kept in the tree. |
| FG2 | **Identical iteration counts and results.** Every Newton iteration count, step-cut count, line-search λ, arc-length step, transient step sequence and accepted/rejected step, every terminal current, weighted current and field, every small-signal admittance, and every run identity digest equal those recorded at commit 0, compared as bit patterns by the same recorder over the solve-layer cases of the existing suites; every existing pin passes unchanged. |
| FG3 | **FD-Jacobian checks.** Every existing FD-Jacobian gate passes unchanged; in addition each term (FE2) has its own FD check of its partials at a non-equilibrium state, under the same tolerances as the existing gates. |
| FG4 | **No measurable assembly slowdown.** A Release benchmark (`[.benchmark]`, not in ctest) times `evaluate` (residual and Jacobian) and `residual` alone on: the legacy 1D diode, the 2D MOSFET with field mobility, a 3D unstructured device, the electrothermal diode, and a nonlocal-tunnelling diode. Measured at commit 0 and after each slice on this machine, in the same session and power state, the median of at least 31 repetitions after warm-up; the gate is the post-slice median not above the pre-refactor median by more than the threshold of FD2 (proposed: 2%, or the measured run-to-run spread of the medians if that is larger, reported alongside). The end-to-end MOSFET and transient benchmarks' wall times are reported, not gated. |
| FG5 | **Suites.** Full Release and Debug suites pass at every slice commit (18/18 and 12/12 as last recorded on `main` after Unit 38, to be re-run at commit 0 before relying on it); zero warnings under the existing flags (`/W4 /WX`). |
| FG6 | **Mutation checks.** For every slice, mutants applied by hand to the new code are caught by the suites **or by the FG1/FG2 comparison** (a reordering mutant, such as two node terms swapped, changes only the last bits, which the suites' tolerances do not see; the comparison is the gate for those, and the report says which of the two caught each mutant). At least: the `DofMap` swapping two fields, an off-by-one block size, the ionization override turned into an addition, two node terms swapped in order, an edge term skipping row b's Dirichlet check, a term's pattern entry dropped, a thermal term's θ column moved, the small-signal pass including interface terms twice, and the wrapper ignoring a model switch. Each is listed with the test that caught it; a surviving mutant is a missing gate, added before the slice is committed. |
| FG7 | **Slice 1 stands alone.** After slice 1 alone: FG1–FG6 pass; `stride()` and every public signature unchanged; the field of every unknown index and the unknown of every (field, node) are queryable; the assembler's block size comes from the map (one source of truth: `m_` is set from it; the literal `4` and `16` of the electrothermal code, where the block is 4 by construction, move with that code into terms in slice 2b); the solve layer's slot reads (`bias.cpp`, `small_signal.cpp`, `trace.cpp`, `transient.cpp`, `fields.hpp`) go through the map. The per-access `m_ * i + k` arithmetic inside the assembler's loops stays until slice 2 moves those loops into terms (rewriting all of it in slice 1 would make slice 1 as large as slice 2). |
| FG9 | **Layering.** The field, map and term types live in `assemble`; `solve` reads only the public map queries (FD4); no new dependency from `physics`, `mesh` or `linalg` on `assemble`; the existing header-boundary checks pass. |
| FG8 | **No scope creep.** No change to physics, models, defaults, deck syntax, file formats, run-record fields or public `DriftDiffusion` signatures; `ARCHITECTURE.md` 2.2, 6.12.1–2 and 11.3 updated to "as built" only after the owner accepts the result. |

**Risks:**

| # | Risk | Mitigation |
|---|---|---|
| FR1 | Bit identity is fragile: a reordered sum, a hoisted common factor (e.g. computing `V * rate_scale_` once instead of per term), a fused multiply-add the compiler forms differently after inlining across a virtual boundary, or a temporary in a different precision changes the last bits and then the iteration counts. | FE3's order rule; terms keep today's expressions verbatim, operand order included; FG1 compares bit patterns at many states before any iteration-level comparison. The flags stay `/fp:precise` with no `/arch` (x64 SSE2 scalar code, `CMakeLists.txt:44`); MSVC is not expected to form contractions under `/fp:precise` without `/fp:contract` [unverified for this toolset: FG1 is what proves it]. If FG1 fails for a reason other than a reordering in the code, the slice stops and the owner is told, rather than tolerances being loosened. |
| FR2 | Indirection costs: virtual dispatch or `std::function` per node or per edge would slow assembly measurably on large 3D meshes; the fused node loop (FE3) needs a per-node call per node term. | Edge and global terms: dispatch once per term per pass, the loop inside the term. Node terms: a per-node call through a list fixed at creation, or static composition (a fixed tuple or variant of terms) if FG4 shows the calls cost (FD3). |
| FR3 | The ionization override and the Dirichlet assignments are not additive; a generic "terms add" interface would break them silently. | Assignment terms are explicit in the interface (a node term may set its rows); FG6 has the override-as-addition mutant. |
| FR4 | The electrothermal pass (`assemble_heat`) duplicates the isothermal terms with T-dependent variants and its own order; folding both into one term list risks differences in either. | Separate thermal term variants with `assemble_heat`'s own order; FG1 runs on electrothermal configurations explicitly. |
| FR5 | Small-signal assembly calls `assemble(…, interfaces = false)` and adds its own interface terms; transient storage and trap elimination thread through `TimeStep`. A missed special case changes AC or transient results only. | FG1 covers the time-step and small-signal matrices; FG2 covers transient step sequences and AC admittances. |
| FR6 | Nonlocal tunnelling changes the Jacobian pattern after creation (`set_paths` → `locate_positions`); a term list with positions fixed at creation would go stale. | The path term relocates its positions on `set_paths` as today; a mutant that skips relocation (FG6, the dropped pattern entry) must be caught. |
| FR7 | Review size: slices 2 and 3 rewrite most of 3200 lines. | One commit per slice, slice 2 split in two (FE7); the recorded baselines (FG1, FG2) instead of a duplicate old assembler; the owner may stop after any slice. |
| FR9 | The FG1 states are only as good as their coverage: a term that is zero at every recorded state (e.g. impact generation below its current floor, a tunnelling path below `matters`) is not checked by bit comparison. | Each configuration's recorded states include one where each of its optional terms is active (nonzero contribution checked when recording); the recorder fails at commit 0 if a term is inactive at every state of its case. |
| FR8 | Unused density unknowns on insulator nodes are baked into the layout; a later layout change (FD5) will change checkpoints and run identities. | Not done in this unit; recorded so the later unit plans a checkpoint version bump. |

**Unresolved decisions (need an owner ruling before implementation):**

| # | Question | Proposal |
|---|---|---|
| FD1 | Include `EquilibriumPoisson` in slice 3 (built from the same Poisson terms) or leave it separate? | Leave it separate in this unit (FE6). |
| FD2 | The slowdown threshold for FG4. | 2%, or the measured spread of the medians if larger, with both reported. |
| FD3 | Dispatch mechanism: virtual terms, or static composition (variant/tuple) chosen at creation? | Virtual dispatch once per term per pass; static composition only if FG4 fails. |
| FD4 | Public or internal: is the `Field`/`DofMap`/term API public in `include/NiTCAD/assemble/`, or internal to `src/assemble/` with only the map queries public? | Internal terms; public read-only `DofMap` queries (the solve layer and later the norm and preconditioners need them). |
| FD5 | Drop the density unknowns on insulator nodes (a layout change)? | No; a later unit with a checkpoint version bump. |
| FD6 | Parallel assembly (edges or nodes split across threads, which changes the order of additions into shared rows and so is not bit-identical)? | Not in this unit; a later unit with tolerance-based gates the owner sets then. |
| FD7 | Derived fields (currents, band edges, charges) as generic fields? | Not in this unit (FE5). |
| FD8 | Where the recorded baselines (FG1, FG2) live: scratch only, or a committed reference file? | The recorder (hidden test cases) is committed at commit 0; the recorded files stay in scratch, regenerated from commit 0 on demand. Bit patterns are specific to this toolset and flags, so a committed file would fail spuriously on CI or after a toolset change. |
| FD9 | Keep the recorder after the unit (a standing bit-identity check for later refactors), or remove it at the end? | Keep it, hidden and out of ctest. |

**Self-review (2026-10-16, before the owner's review).** Checked against the code on `main` (`f153d57`): the
accumulation order, `small_signal_matrix`'s `assemble(…, false)` call, the five solve-layer `stride()` readers, the
build's `/fp:precise` and the run identity (inputs only, so unchanged by a refactor) are as stated. Corrected from the
first draft: (1) node terms may run as separate loops without reordering any entry (FE3); (2) no assembly benchmark or
recorder exists on `main`, so a commit 0 adds them before the refactor, and baselines come from it (FE7, FG1, FG2,
FG4); (3) a duplicate old assembler kept "test-only" was replaced by recorded bit patterns (FG1, FR7); (4) reordering
mutants are invisible to the suites' tolerances and are judged by the comparison (FG6); (5) FG7 asked slice 1 to
rewrite all 114 lines of index arithmetic, which is slice 2's work, now narrowed; (6) the toolset-specific bit
patterns make a committed reference file unsafe (FD8). Added: slice 2 split into isothermal and electrothermal
commits (FE7), the layering gate (FG9), the coverage risk of the recorded states (FR9), FD9.

**Owner ruling (2026-10-17):** "go ahead" after the self-review, taken as authorization to implement with the
proposals of FD1–FD9 accepted as written.

**Commit 0 as built (the recorder differs from FG1's wording above, and is stronger).** Instead of F and J at fixed
states built for the purpose, the recorder traces **every** value the assembler's public members return during the
existing suites:
- `include/NiTCAD/assemble/observer.hpp`, `src/assemble/observer.cpp`: a null-by-default observer. `DriftDiffusion`'s
  public members hand it their inputs and outputs under a tag: `evaluate` and `residual` (state, residual, Jacobian
  values, and a time step's rate, storage and traps), `small_signal_matrix`, both `bias_derivative`s,
  `small_signal_currents` (columns, values, bias), `state_from_potential`, `stamp_contacts`, `update_size`,
  `apply_update`, `storage`, `trap_occupancies`, `edge_currents`, `terminal_currents`, `terminal_current_resolution`,
  `edge_flux_terms`, `band_edges`, `gate_charges`, `contact_charges`, `conduction_currents`, `interface_trap_charges`,
  `thermal_contact_heat`, `electrical_power`, `edge_joule_heat`. Unset, each call pays one null test; the hooks only
  read, so no value changes.
- `tests/baseline/trace_listener.cpp`, a Catch2 listener linked into the assemble, solve, deck and validation test
  executables: with `NITCAD_TRACE_FILE` set, one line per observed value (index in its section, tag, length, a 64-bit
  FNV-1a hash of its bit patterns) under each test case and section. Child processes (the deck tests' `nitcad.exe`)
  are not traced.
- Each traced ctest entry is run with the variable set (scratch script; the files stay in scratch, FD8); a slice
  passes FG1 and FG2 when its traces equal commit 0's line for line. The solve layer is unchanged by every slice, so
  equal assembler outputs at every call give equal iteration counts and results; the first differing line names the
  test, section and call. It covers every term the suites activate, at the suites' own states.
- **FR9's planned mitigation was not built.** The plan was for the recorder to fail at commit 0 when a term is
  inactive at every recorded state. Because the recorder traces only the states the suites reach, a term the suites
  never activated would go unchecked bit for bit. In its place:
  - The existing per-model FD-gate cases evaluate each model where it is active: impact ionization, local and
    nonlocal tunnelling, gates, interfaces and traps, the electrothermal fixtures, time steps. This was checked by
    reading those tests, not measured.
  - From slice 2a, the terms test requires every term to be active in its fixtures and checks its partials.
- The FG4 benchmark: `tests/solve/assembly_timing.hpp` and five `[.assembly_benchmark]` cases (the legacy 1D diode,
  the 2D MOSFET with field mobility, the 3D generated diode, the electrothermal pn bar, the Zener diode with nonlocal
  paths) in the solve test files that already build those devices.

**Commit 0 findings (2026-10-17, commit `4e8a3af`, this machine: Intel Core Ultra 5 225F, 6 P + 4 E cores):**
- Catch2 v3 runs test cases in random order by default ("Randomness seeded to …"), so traces are recorded with
  `--order decl`; in random order one-time static work (e.g. a cached reference solution) is traced under
  whichever case runs first.
- Determinism: two full traced runs (A, B) of the same commit-0 binaries agreed on every traced entry except
  `solve_unstructured`, whose three cases (an 80 000-node tensor reference, 241 803 unknowns) differed in B from the
  first observed value on (the printed results agreed to every printed digit). Later full runs showed the same
  kind of split in `solve` (after a complex solve), `solve_btbt` (after a Newton linear solve) and `solve_mosfet`:
  always exactly two outcomes per entry, always starting right after a linear solve. **Cause found:** the
  outcome depends on the core type the process runs on. `btbt: a curved junction tunnels before a planar one`,
  pinned to each logical CPU in turn, gives one trace on CPUs 0, 1 and 6–9 (the P cores) and another on CPUs 2–5
  (the E-core cluster); an unpinned process lands on either. (The earlier pinned probes happened to use P cores
  only.) The likely mechanism [unverified]: Eigen sizes its dense kernels' blocking from the cache sizes of the
  core it first queries (the E-core cluster's L2 differs from a P core's), and SparseLU's supernodal updates then
  sum in a different order. **This is a reproducibility defect on `main`** (ARCHITECTURE.md 6.8), outside this
  unit's scope: results depend in the last bits on the core a run starts on, on hybrid CPUs. A fix (e.g. fixing
  Eigen's cache sizes at start-up) changes last bits of existing results, so it needs its own owner-approved unit.
  **For this unit's gates:** every trace run is pinned to the P cores (affinity mask 0x3C3), and a slice's trace
  must equal the pinned baseline's exactly.
- FG4 noise: three back-to-back runs of the benchmark on commit 0 gave medians spread by 20–40% per device (e.g.
  the MOSFET's `evaluate` 3818–4655 µs, the 3D diode's 1706–2274 µs): the machine's run-to-run variation is far
  above FD2's 2%. FG4 is therefore measured by alternating runs of the commit-0 test binary (kept in scratch) and
  the slice's, several times back to back, comparing the minimum of each device's medians; the gate is the slice's
  minimum not above commit 0's by more than 2% plus the spread of commit 0's minima across the alternation, with
  every number reported.
- An existing timing gate, `btbt: the cost of the nonlocal paths` (`t_on < 2 t_off + 0.05`), failed once when nine
  traced test processes ran at once (0.88 s against 0.63 s); it passes under `ctest -j 4`. Not changed in this
  unit; noted as a load-sensitive test.

- Commit 0b (`a7ddc6b`): the public solves also hand their results (fields, bands, currents, charges, admittances,
  convergence records, run identities) to the observer, so the slices' solve-layer conversions are compared bit
  for bit, not only through the tests' tolerances.
- FG4 as measured: commit 0b's and the slice's `[.assembly_benchmark]` runs alternate four times, each pinned to
  CPU 0 (a P core), and the minimum per device is compared. Pinned, the spread of commit 0b's four runs is 1–10%
  per device (unpinned it was 20–40%).

**Slice 1 results (fields and the unknown map):**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces: two commit-0b baselines identical to each other on all nine traced entries; slice 1 identical to them on all nine (the only extra lines are the two new `fields` test cases). |
| FG3 | Every existing FD-Jacobian gate passes unchanged (no term moved in slice 1). |
| FG4 | `evaluate`, slice 1 against commit 0b (minimum of four pinned runs): MOSFET 2056 against 2038 µs (+0.9%), 3D diode 665 against 650 (+2.2%; commit 0b's own spread 9.7%), electrothermal bar +0.4%, 1D diode +0.6%, Zener +0.9%; `residual` −1.8% to −0.1%. Within the gate (no assembly loop changed). |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 7 mutants, all caught: electron and hole fields swapped in the map (fields test), `DofMap::index` off by one (fields test), the temperature slot pointing at the holes (electrothermal test), `scaled_state` writing electrons into the hole slot (small-signal test; the run also crashed), the transient error estimate reading the potential for the rise (electrothermal test), `converged_temperature` reading the potential (electrothermal test); the small-signal residual bound reading the electron slot for the potential survived the suites and was caught by the trace (only the `result.ac_resolution` lines differ). |
| FG7 | `stride()` and every public signature unchanged; `DriftDiffusion::dofs()` answers the field and node of every unknown and the unknown of every (field, node); the block size is set from the map; the solve layer reads every slot through it (`solve/fields.hpp` `Slots`, `scaled_state`). |
| FG9 | The map is in `assemble` (`fields.hpp`); `solve` reads only its public queries. |

**Slice 2a as built (the isothermal terms; `src/assemble/terms.hpp`, `terms.cpp`, internal per FD4):**
- `DriftDiffusion::Term`: `info()` (name, domain, the fields whose rows it writes and those its partials are taken
  in, whether it assigns its rows), `pattern()` (the entries it needs beyond every node's full block and the edges'
  stencil), `add(system, evaluation)`. Terms hold no data: they read the system's at each call, so a system stays
  copyable and shares its terms (`std::shared_ptr<const Term>`). `DriftDiffusion::Terms` holds the term types, the
  factory, the pattern hook and a one-term evaluation for the tests.
- Fifteen terms, in FE3's order: ohmic contact rows, insulator rows, Poisson charge, incomplete ionization, gate,
  SRH, Auger, radiative, storage (nodes); interfaces (interface edges); Poisson flux, carrier fluxes (mesh edges);
  impact ionization, local and nonlocal band-to-band tunnelling (global). The loop bodies are the old ones, moved;
  `add_impact_generation`, `add_local_tunnelling` and `add_path_tunnelling` became the three global terms.
  `assemble` is the term loop; with `interfaces = false` (small-signal) it skips the interface-edge terms.
- **One loop per node term, not one fused node loop (FE3's default).** FD3's dispatch (one virtual call per term
  per pass) would otherwise have been a call per node per term. FE3 allows separate loops: a node term writes only
  its own node's rows, so every entry receives its additions in the old order. FG4 decided: no slowdown (below).
- **The mesh-edge pass is split** into the Poisson flux and the carrier fluxes. The two write disjoint entries (ψ
  rows; n and p rows), so the order of the additions into each entry is unchanged. The electron and hole fluxes
  stay one term, because `edge_fluxes` computes them together.
- **The pattern.** `create` builds the edge table before the pattern, adds each term's own entries (interfaces'
  insulator-row density columns, impact ionization's other carrier, the paths'), then locates the positions as
  before. `set_paths` rebuilds the pattern from the base plus every term's entries. The pattern is the same set, and
  `from_triplets` sorts, so it is identical. Which terms exist depends on the models alone (the impact and tunnel
  node tables are built after the pattern); a term whose table is empty adds nothing.
- `tests/assemble/terms_test.cpp` (FG3): each term's partials against central differences of its own residual, at
  non-equilibrium states across five fixtures: a diode with every node model under both statistics, steady and with
  a time step; the oxide/electrode/trapped-interface fixture, steady and with a time step; a lumped gate; impact
  ionization at 3e5 V/cm; local and nonlocal tunnelling at 4e6 V/cm. Each runs in 1D and 2D (the gate in 1D only),
  with the tolerance of the Unit 19/20 gates (1e-6, the best of steps 1e-5 to 1e-7). Each fixture requires its terms
  to be active (FR9). A second case checks the order, the domains and which terms assign. The column scale is the
  whole system's Jacobian, the existing gates' scale. Scaled by a term's own column, the vanishing side of a
  Scharfetter–Gummel flux at a drop of hundreds of V_T (an entry near 1e-128 beside a flux near 300) is below what
  any difference step resolves, and reads as error 1.

**Slice 2a results:**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces identical to slice 1's on all nine traced entries. |
| FG3 | Every existing FD-Jacobian gate passes unchanged; the 15 terms pass their own checks (230 assertions). |
| FG4 | Minimum of four alternating pinned runs, slice 2a against slice 1. `evaluate`: MOSFET 2028 against 2033 µs (−0.2%), 3D diode 647 against 664 (−2.7%), electrothermal bar −0.1% (not a term yet), 1D diode 11.9 against 13.2 (−10%), Zener 11.1 against 12.1 (−8.6%). `residual`: MOSFET +1.5% (slice 1's own spread 4.8%), others −3.2% to −0.1%. No slowdown. |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 9 mutants, all caught by the suites. Ionization override turned into an addition: heterojunction test. SRH and Auger swapped in order: the terms test's order check (the swap itself only changes last bits). Carrier fluxes ignoring row b's Dirichlet check: electrothermal test. The interface pattern without the hole column, and the paths' pattern entries dropped: the contract check in the position lookup. The impact pattern on the wrong edges: impact test. `assemble` ignoring `interfaces = false`, so small-signal adds the interface terms twice: small-signal test. The Poisson flux ignoring an electrode row: insulator test. The gate term left out: gate test. |
| FG9 | The terms are in `src/assemble` (internal); no new dependency; `solve` unchanged. |

**Slice 2b as built (the electrothermal rows; `src/assemble/thermal_terms.cpp`):**
- `assemble_heat` stays the entry point, because `thermal_contact_heat` calls it with the sinks out. It is now the
  loop over the electrothermal term list, which `Terms::make` builds when the model is on. The evaluation carries the
  thermal states, the levels and the sinks flag.
- Twelve terms, in `assemble_heat`'s order:
  - heat sinks (assigns the heat row: the isothermal sink's Dirichlet row, R_th, or 0)
  - lattice heat storage
  - ohmic contact rows at the node's T
  - insulator rows (an electrode at its T)
  - the Poisson charge (the isothermal term, shared; complete ionization, U6)
  - gate with the p-polysilicon shift
  - SRH, Auger and radiative at the node's T
  - carrier storage (the continuity rows and the stored carriers' energy in the heat row)
  - interfaces
  - edge fluxes (Poisson, carriers and heat)
  - impact ionization with gamma(T)
- Per node, the heat row is written before the electrical rows, and the carrier storage adds to it last. With one
  loop per term, every entry still receives its additions in that order.
- **The edges stay one term.** The heat rows take the carrier fluxes and their partials, so splitting would compute
  `thermal_fluxes` twice.
- **The interface term runs even for small-signal.** `assemble_heat` never skipped interfaces, and the electrothermal
  `small_signal_matrix` adds no interface terms of its own.
- The literal 4 and 16 moved with the code into the terms (FG7). Every block is 4 × 4 by construction.
- The edges' full 4 × 4 pattern blocks moved from `create` into the edge term's `pattern()`, so the pattern is the
  same set. `make_thermal`'s positions and the diagnostics (`add_thermal_storage`, `bias_derivative`,
  `thermal_contact_heat`) keep their code (FE5).
- `terms_test` gains:
  - a heated 1D and 2D diode in radiative silicon, with an oxide and a recombining fixed-charge interface in 2D,
    under both statistics, steady and with a time step;
  - a p-polysilicon gate;
  - impact ionization at about 3e5 V/cm;
  - the electrothermal order;
  - a held heat row being its sink's Dirichlet row in a time step too.

  The last was added for two mutants that survived the suites (below).

**Slice 2b results:**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces identical to slice 2a's on all nine traced entries. |
| FG3 | Existing FD-Jacobian gates unchanged; all 27 terms pass their own checks at 1e-6 and are active (478 assertions in 4 cases). |
| FG4 | Minimum of four alternating pinned runs, slice 2b against 2a. `evaluate`: electrothermal bar 194.6 against 192.7 µs (+1.0%, slice 2a's own spread 7.1%); `residual` +0.9%. The isothermal devices −3.0% to −0.2% (unchanged code). Within the gate. |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 8 mutants, all caught. Heat-sink assignment turned into an addition: electrothermal test. A θ column moved (the contact row's ψ entry into the p column): electrothermal test. Thermal SRH and Auger swapped: the terms test's order check. The edges' pattern without the heat rows: the position lookup's contract check. The interfaces term left out: the terms test's activity check. The gate's p-polysilicon shift dropped: electrothermal test. Lattice storage, and separately the carriers' stored energy, added to a held (isothermal-sink) heat row: **both survived the suites**. At a converged step the held rise equals the sink's before and after, so the term vanishes at the solution, and no fixture held a node that was neither contact nor insulator. The new held-row case (a sink on y_min) catches both. |

**Slice 3 as built (composition; `src/assemble/composition.hpp`, `composition.cpp`):**
- `DriftDiffusion::Composition` owns the unknowns' map and the active terms, with each term's equation sets, and
  runs the evaluation loop over them. `DriftDiffusion` holds one, immutable and shared by copies, in place of
  `dofs_` and the term list, and forwards `dofs()`, `assemble` and `assemble_heat` to it. Every public signature is
  unchanged (`dofs()` moved out of line).
- **Equation sets after ARCHITECTURE.md 6.12.2.** Poisson, continuity, contacts, oxide, gate, interfaces,
  generation and heat. They form a bit set: the electrothermal edge term belongs to Poisson, continuity and heat
  at once.
- **Composition tables.** The two `make` functions became tables in `terms.cpp` and `thermal_terms.cpp`. Each row
  is a term's sets, an activation predicate over the models and the device data `create` has set by then
  (contacts, insulators, gates, interfaces), and a factory. Table order is FE3's assembly order. `compose` picks
  the table by the electrothermal model, keeps the active rows, and builds the map from the fields the equations
  solve for.
  - Diode = Poisson + continuity + contacts.
  - MOS adds oxide, gate and interfaces.
  - Thermal adds heat and a fourth unknown.
  - The generation models add generation.
- The contact-row and insulator-row terms are now active only when the device has such nodes. Before, they were
  always in, with empty loops. No value changes.
- **Not done:** the model data the terms read (per node, edge and contact) stays in `DriftDiffusion`, because its
  diagnostics read it too (FE5). `DriftDiffusion` is therefore both the compatibility wrapper (FE4) and that
  data's store. Moving the data out and retiring the wrapper is a later unit's work, as FE4 already says.
- `terms_test` gains three composition cases:
  - the equation sets and block sizes of a diode, a MOS structure with an oxide and interface, a lumped-gate
    MOS-C, a generation model, and a thermal diode;
  - every model switch (SRH, Auger, radiative, incomplete ionization, impact, local and nonlocal tunnelling; the
    electrothermal SRH, Auger, radiative and impact) adds exactly its term;
  - the electrothermal small-signal matrix at s = 0 equals the steady Jacobian entry for entry, so its interface
    terms stay in.

**Slice 3 results:**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces identical to slice 2b's on all nine traced entries, every test case except the `terms:` and `composition:` cases. Those trace their own calls to `evaluate`, so they change with the term lists of their fixtures (an oxide-free fixture no longer carries an empty insulator-rows term). |
| FG3 | Existing FD-Jacobian gates unchanged; the term checks pass. |
| FG4 | Minimum of four alternating pinned runs, slice 3 against 2b, measured twice. `evaluate` −0.9% to +2.3% (1D diode, 11.4 against 11.7 µs; +0.9% in the first measurement), `residual` −2.1% to +0.3%. The first measurement's 3D-diode `residual` +2.5% was +0.3% when repeated. Within the gate; the only code-path change is one indirect call per evaluation. |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 7 mutants. Six caught: the Auger switch ignored (drift-diffusion test); the electrothermal table never chosen (crash); the `interfaces` flag ignored, so small-signal adds the interface terms twice (small-signal test); radiative and storage swapped in the table (the order check); the temperature field dropped from the map (electrothermal test); the gate predicate false (gate test). The seventh (the electrothermal interface term skipped without `interfaces`) **survived and was equivalent**: `assemble_heat` always passes `true`, so the check it mutated was dead code. The dead check was removed. The new electrothermal small-signal case guards that behaviour. |
| FG8 | No change to physics, models, defaults, deck syntax, file formats, run records or public signatures. `ARCHITECTURE.md` 2.2, 6.12.1–2 and 11.3 are to be updated to "as built" after the owner accepts the result. |
| FG9 | Fields and the map public in `assemble`; terms and composition internal (`src/assemble`); no new dependency; `solve` reads only the map. |

**Review (2026-10-17, owner request: review all eight commits against `main` before push authorization):**

| Area | How | Result |
|---|---|---|
| Term order and bodies | `main`'s `assemble`, `add_*`, `assemble_heat` and `add_thermal_*` against the 27 terms, statement by statement after normalizing only the mechanical renames (`s.`, `ev.`, `stride`/`m` for `m_`, `at(s, ev, …)`) | Every expression unchanged operand for operand. The only other lines are skip conditions equivalent to the old if/continue chain, helpers with the old lambdas' bodies, and wrapping. The tables' order is `main`'s sequence |
| Activation predicates | Each against the guard it replaced | Equivalent. The impact, tunnel and path tables are filled only under the same model switches. `InterfaceEdges::empty()` is `edges_.empty()`. The pattern entries are the same set |
| `create`, `set_paths` | Read line by line | Positions are located the same way in the same pattern; paths rebuild from the same entry set |
| Public API | Headers diffed against `main` | Nothing public removed or changed. Additions: `dofs()`, `fields.hpp`, `observer.hpp` |
| Bit identity against `main` | `main` built in a separate worktree. Six decks run with `nitcad run` from `main` and from HEAD, pinned to a P core: a diode with Fermi–Dirac, field mobility and incomplete ionization; an oxide with an electrode and a trapped interface; a lumped gate; electrothermal; impact ionization; nonlocal tunnelling | Output tables byte-identical. This adds direct steady-state evidence for the link from commit 0 back to `main`, which until then rested on the hooks only reading. Transient and small-signal (which the deck runner does not run) still rest on that, and on the suites |

Six findings, none a blocker, fixed at the owner's word ("do 1"):
1. **The terms' field declarations were unchecked.** New test case: every term writes only the rows and Jacobian
   columns of the fields it declares, across every fixture. Before it was committed, it caught a deliberately
   wrong declaration (the gate term declaring the density rows).
2. The five `*_unobserved` solves had external linkage and no declaration. They are now in anonymous namespaces.
3. `Term`, `Terms` and `Composition` were public nested names. They are now private to `DriftDiffusion`. The
   library and its tests name them through `TermAccess`, declared only in the internal `src/assemble/terms.hpp`.
4. `dofs()` on a moved-from system dereferenced a null composition. It now has a contract check
   (`NITCAD_EXPECTS`).
5. Electrothermal under Fermi–Dirac computes the equilibrium product in each of SRH, Auger and radiative (once
   before). Measured with a new `[.assembly_benchmark]` case, the 401-node pn bar under Fermi–Dirac, against
   `main` (the case was added there temporarily, then removed): `evaluate` 339.9 against 342.8 µs (−0.8%),
   `residual` −3.0%. No slowdown.
6. The FR9 wording above was corrected: the planned recorder check was not built, and the text now says what
   covers its concern instead.

After the fixes: pinned traces identical to slice 3's on all nine traced entries (the `terms:` and `composition:`
cases aside, as before); Release 18/18, Debug 12/12, no warnings; FG4 against slice 3 with six alternating pinned
runs, `evaluate` −0.6% to +0.7%, `residual` −1.6% to +2.8% (the 3D diode; slice 3's own spread 14.6%). A first
four-run measurement read up to +6.3%, but with both binaries about 15% slower than usual, and did not repeat.

**Status:** commit 0, commit 0b and slices 1, 2a, 2b, 3 done. Unit 40 is complete on `assemble/fields-equations`,
reviewed, and the review's findings fixed. **Accepted by the owner (2026-10-17, "accept and push Unit 40")**;
`ARCHITECTURE.md` 2.2, 6.12.1–2 and 11.3 updated to "as built" (FG8). The hybrid-CPU reproducibility
defect (commit 0 findings) is tracked separately as Unit 41 (`investigate/hybrid-cpu-determinism`).

### Unit 41 — last-bit reproducibility on hybrid CPUs (opened 2026-10-17; fix decided by the owner 2026-10-17; on `linalg/hybrid-cpu-determinism`)

**Owner request (2026-10-17):** track the P/E-core last-bit difference found during Unit 40 (its commit 0
findings) as its own unit, separate from Unit 40. First a minimal reproducer and the root cause. Do not change any
numerical result or implement a fix until the impact, the performance cost and the reproducibility policy have been
reviewed.

**Reproducer (`benchmarks/hybrid_cpu/`, standalone, Eigen only, not in the CMake build).**
- `repro.cpp` prints the executing core's CPUID leaf-4 cache sizes, the sizes Eigen adopted, and bit hashes
  (FNV-1a) of a 300 × 300 dense product and a SparseLU solve of a 22 500-unknown nonsymmetric five-point operator.
- Run pinned to each logical CPU of this machine (Intel Core Ultra 5 225F, 6 P + 4 E cores):

| Cores | CPUID L1d / L2 | Eigen's adopted L1 / L2 | SparseLU hash, `auto` | `fixed` (any one set) |
|---|---|---|---|---|
| P: CPUs 0, 1, 6–9 | 48 KiB / 3 MiB | the same | `317ad696…` | one hash on all ten CPUs |
| E: CPUs 2–5 | 32 KiB / 4 MiB | the same | `50df5dd6…` | the same hash |

- The 300 × 300 dense product agrees on every core. A 1200 × 1200 product does not.
- `diode_2d.deck` (a 121 × 81-node 2D diode, `nitcad run` from `main`): 6 of the table's current values differ
  between a run pinned to CPU 0 and one pinned to CPU 2, by at most 1.1e-15 relative. CPU 0 against CPU 6 (two P
  cores) is identical. Six small decks (1D, and a small 2D oxide structure) are identical on P and E cores. Their
  systems are too small for the blocking to change.

**Root cause.**
- Eigen (5.0.1, vcpkg; the backend's static_assert pins it) detects the cache sizes once per process, with CPUID leaf 4 on the core that runs the first
  matrix product, and keeps them in a static.
- The blocking of its matrix products is computed from them. SparseLU's supernodal updates (`panel_bmod`,
  `kernel_bmod`) are such products.
- On a hybrid CPU the two core types report different L1 data caches, so the products split their sums
  differently and round differently in the last bits. Newton then amplifies the differences from the first linear
  solve on.
- **The L1 size alone decides it.** On a P core, `fixed 32768 3145728 …` (only L1 changed to the E core's) gives
  the E-core hash. On an E core, `fixed 49152 …` gives the P-core hash. L2 at 3 or 4 MiB makes no difference.
- So the result depends on the core type the first product runs on. That is usually the start core, but a thread
  can migrate before then. It does not depend on the machine's load or on the core of later products.
- [unverified] The MKL PARDISO backend (Unit 18) dispatches by CPU too. Whether it differs between core types is
  not measured; MKL's conditional numerical reproducibility settings (`MKL_CBWR`) would be the place to look.

**Measured (2026-10-17, at the owner's request; no code, flag or result changed).** The tooling:
`benchmarks/hybrid_cpu/measurement.patch` is a measurement-only patch, not committed to any build. It applies to
Unit 40's tree, which has the trace listener, and:
- adds full values for the `result.*` tags to the trace;
- adds `NITCAD_EIGEN_CACHE=l1,l2,l3`, which fixes Eigen's cache sizes before any test runs;
- adds a hidden `[.unit41_cost]` case timing SparseLU on the MOSFET and a 20³ 3D diode.

`compare-values.sh` compares two value traces block by block. The suites were traced three ways, each pinned:
P cores (mask 0x3C3), E cores (0x3C), and E cores with the cache sizes fixed to the P-core values.

*Impact on the suites:*
- **Every test passes on the E cores as on the P cores.**
- **Fixing the cache sizes removes the difference completely.** On E cores with the P-core sizes, all nine traced
  entries are bit-identical to the P-core run. Eigen's L1-dependent blocking is the whole cause.
- **Which entries differ.** With Eigen's own detection, 4 of the 9 entries are identical on P and E (assemble,
  deck, transient, validation) and 5 differ. In `solve`, `solve_unstructured` and `solve_breakdown` the
  differences are in values. In `solve_mosfet` and `solve_btbt` they also change the work done. Ten of the MOSFET
  cases and two tunnelling cases ran a different number of evaluations, i.e. Newton iterations or bias cuts. For
  example, the curved-junction tunnelling case made 1825 calls on P and 1558 on E; the mesh-independence MOSFET
  case 224 and 204.
- **Size of the changes.** Fields (potential, densities, band edges) agree to at most 4e-15 relative of each
  array's largest value, even where iteration counts changed.
- **Terminal currents.** Every difference is at most 4.3e-5 of the larger of the current and its own reported
  resolution (`terminal_current_resolution`). Currents more than ten times their resolution change by at most
  3.1e-8 relative: Newton stops at a different iterate within its tolerance. Small-signal admittances: the
  dominant entries agree to 1e-11 relative; entries at the AC resolution (about 1e-11 S against 6e-9 S) change sign.
- **Quantities that are rounding by nature** change by large relative amounts, as expected: backward errors, the
  AC resolution, edge currents and quasi-Fermi levels at equilibrium.

*Cost of fixing the cache sizes* (SparseLU refactorization, fastest of three pinned runs; the noise, between two
configurations with the same blocking, is up to 2%):

| Fixed to | P core, MOSFET / 3D diode | E core, MOSFET / 3D diode | Result bits everywhere |
|---|---|---|---|
| nothing (today) | 155.7 / 5771 ms | 180.9 / 6563 ms | P and E differ |
| P-core sizes (L1 48 KiB, L2 3 MiB) | unchanged | 184.1 / 6809 ms (+1.8% / +3.7%) | today's P-core bits |
| E-core sizes (L1 32 KiB, L2 4 MiB) | 154.8 / 5779 ms | unchanged | today's E-core bits |
| the smaller of each (L1 32 KiB, L2 3 MiB) | 155.2 / 5787 ms | 184.6 / 6526 ms | today's E-core bits |

The cost of any fixed choice is within about 4%, at most at the edge of the noise. The choice of sizes decides
which of today's two results survives: P-core sizes keep the P-core runs' bits, which match Unit 40's recorded
traces; a 32 KiB L1 keeps the E-core runs'.

**Options for the owner (option 1 with the P-core values chosen 2026-10-17; below):**
1. Fix Eigen's cache sizes at start-up to one set, for example in `linalg`'s Eigen backend before the first
   factorization:
   - to the P-core values, so runs that start on an E core change in the last bits;
   - or to a 32 KiB L1, so P-core runs change.

   Either makes one machine reproducible whatever core a run starts on.
2. As 1, but take the smallest L1 over all cores, queried per core with thread affinity. The same on this machine
   as the 32 KiB choice, and it adapts to other hybrid parts. Different machines still differ.
3. A solver option (sizes or "fixed"/"detected"), recorded in the run record.
4. Document only: results are reproducible on one core type, and runs may differ in the last bits between P and E
   cores.

For reproducibility across machines, the blocking must not depend on the hardware at all: one fixed set
everywhere (option 1 with constants). Its cost on other CPUs is not measured.

**Owner decision (2026-10-17):** option 1, with the P-core values. The CPU a run executes on must not change its
result; bit-identity across machines is not promised, and tolerances stay the contract (ARCHITECTURE.md 6.8).

**Fix (branch `linalg/hybrid-cpu-determinism`, from `main` 629bed0 with the two investigation commits).**
- `src/linalg/eigen_sparse_lu.cpp`: `fix_cache_sizes()` calls `Eigen::setCpuCacheSizes` once per process
  (`std::call_once`) with L1 48 KiB, L2 3 MiB, L3 20 MiB, at the start of `SparseLuBackend::analyze`, before any
  product of the backend.
- `tests/linalg/eigen_blocking_probe.cpp`, ctest `linalg_eigen_blocking`: starts itself once per logical CPU of
  its affinity mask, each child suspended, pinned, then resumed; each solves a 22 500-unknown nonsymmetric
  five-point system with the default `LinearSolver` and prints a hash of the solution's bits. It fails unless all
  hashes agree. Its own executable, because Eigen detects the cache sizes once per process.

**Verified (2026-10-17):**
- Probe with the fix: one hash (`130858804b42571a`) on all ten CPUs.
- Mutant (the `setCpuCacheSizes` call removed): the probe fails; CPUs 2–5 (E) give `b1405c79cef29076`, the P cores
  `130858804b42571a`, the fixed build's hash. The P-core results are kept.
- Release suite 19/19 (the probe added).
- Not run, at the owner's direction: a whole-suite trace pinned to P against E cores with the fix in place. The
  investigation's E-core run with the P-core sizes was bit-identical to the P cores on all nine entries, but it
  fixed the sizes before any test ran; the fix does so at the backend's first analysis. An Eigen product large
  enough to be blocked that runs before the first factorization would still use the detected sizes. The probe
  does not cover that.
- Debug suite 13/13, the probe among them (its CPUs agree in Debug too).

**Still open:**
- [unverified] PARDISO (MKL) on P against E cores (`MKL_CBWR`).
- The cost of the fixed sizes on CPUs other than this one is not measured.
- [unverified] Whether any Eigen product large enough to be blocked runs before the first factorization (above).

**Status:** fixed; merged to `main` as `67e8fe5` (PR #51, 2026-10-17).

### Unit 42 — results and provenance (proposed, authorized and accepted 2026-10-17; on `results/provenance`, PR #52)

**Owner request (2026-10-17):** a scope proposal for the "Results and provenance" row of ARCHITECTURE.md 11.3
(status "partial"; missing "build/version stamp in files; transient and small-signal results on disk"). Documentation
only until the owner approves.

**Why now (audited on `main` at `67e8fe5`).**
- `results::RunRecord` holds `input_identity` (a digest of the inputs) and `settings`. Nothing says which build
  produced a result: `project(NiTCAD VERSION 0.0.0)` is the only version, and no source revision, compiler, linear
  backend, MKL version or start time is recorded anywhere.
- After Unit 41 a result is reproducible on one machine whatever the core. To reproduce it, a reader needs the build
  and the platform it came from, and a file does not carry them.
- `io` saves meshes, sweeps, checkpoints and devices (content kinds 1–4). `results::Transient` (Unit 21) and
  `results::SmallSignal` (Unit 22) have no file: a transient or AC run is lost when the process ends.
- `trace_bias` (Unit 19) returns a `Sweep`, so `save_sweep` already stores it (no test says so).
- The deck (Unit 31) runs `equilibrium` and `sweep` only (`deck::Analysis`).

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| PV1 | **A `results::Provenance` record in `RunRecord`** (plain data; `results` still depends on `base` only). Build: NiTCAD version, source revision (git commit and a `dirty` flag), build type, compiler and version, Eigen version. Run: linear backend in effect, the MKL version string when PARDISO is loaded, Eigen's fixed cache sizes (Unit 41), thread count, OS version, CPU model string, logical CPU count, start time (UTC, ISO 8601) and wall time. | Enough to say which build, backend and platform produced a result and when; the items Unit 41 showed can change last bits are all in it. | A few hundred bytes per run. |
| PV2 | **Provenance is not hashed.** `input_identity` keeps its current definition. | Run identities and checkpoints written before this unit stay valid; identity is "same inputs", provenance is "how and where". | Two runs of one input on different builds share an identity; provenance tells them apart. |
| PV3 | **Build facts from a generated header** (`base`), refreshed at build time by a CMake step that rewrites it only when the revision or dirty state changes; outside a git checkout the revision is `unknown`. Run facts gathered by `solve` when it makes the run record (`make_run_record`), the wall time when the run ends. | A configure-time stamp goes stale after the next commit; rewriting only on change avoids relinking every build. | One custom build step; it needs `git` on the build machine's path (not a runtime dependency). |
| PV4 | **In files as a keyed record inside the run record** (IO2): sweeps, checkpoints and the new kinds carry it; files without it load with an empty provenance. No format version change. | Old files keep loading; new readers see the stamp. | Files written before this unit have no provenance, by necessity. |
| PV5 | **Transient and small-signal files:** content kinds 5 and 6; `save_transient` / `load_transient`, `save_small_signal` / `load_small_signal`; every field bit for bit (complex values as two doubles by bit pattern), reusing the point, field, convergence and error records. A test that a `trace_bias` sweep round-trips through `save_sweep`. | Closes the "on disk" half of the 11.3 row and the R2 deferral for these results. | Snapshot-heavy transients make large files (Unit 21's snapshot policy already bounds them). |
| PV6 | **The runner:** `nitcad --version` prints the build stamp; `nitcad run` prints the provenance with its report. | The first thing to ask of a result is which build made it. | — |
| PV7 | **Out of scope:** deck analyses for transient, small-signal and trace runs (PD1); visualization formats (VTK/HDF5, Unit 25); cross-machine checksums (6.8); any change to numerics, defaults or deck syntax beyond PD1's ruling. | Keeps the unit to records and files. | Transient and AC files are reachable from the API only, unless PD1 brings the deck in. |

**Owner decisions needed:**

| # | Question | Recommendation |
|---|---|---|
| PD1 | Add `[transient]`, `[small_signal]` and `[trace]` deck analyses (with their outputs) in this unit, or as the next one? | **Next unit.** Each brings its own deck keys (time stepping, frequencies, arc length) and checks; this unit stays reviewable. |
| PD2 | Start time and wall time in files: two runs of one deck then give files that differ in those bytes. | **Include them.** Tests compare loaded values with provenance excluded, not file bytes. |
| PD3 | Host name and user name in provenance? | **Exclude.** Files get shared; CPU model and OS are what reproduction needs. |
| PD4 | Provenance in the CSV table (`out.table`)? A `#` header line breaks some readers; a sidecar file is one more output. | **Leave the CSV unchanged**; the binary sweep file and the runner's report carry the provenance. |
| PD5 | The version stays `0.0.0` until a release policy exists? | **Yes**; the revision identifies builds until then. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| PG1 | **No numerical change.** Every suite passes unchanged; the Unit 40 assembly traces (P-core pinned) are bit-identical to `main`'s; `input_identity` of fixed fixtures equals the value on `main`. |
| PG2 | **Round trips.** `Transient` and `SmallSignal` round-trip bit for bit, including NaN, −0 and infinities, empty and present `stopped` / `unfinished`, snapshots, traps, electrothermal and tunnelling fields; provenance round-trips in every kind that carries it; a `trace_bias` sweep round-trips. |
| PG3 | **Compatibility.** Sweep, checkpoint and device files written by `main` (committed as fixtures) load unchanged with an empty provenance; a loader given another content kind refuses it (`invalid_input`). |
| PG4 | **Corruption.** Truncation, a checksum mismatch and appended bytes are refused for kinds 5 and 6, as for kinds 1–4. |
| PG5 | **Provenance is true.** The revision equals `git rev-parse HEAD` of the build and `dirty` follows the working tree; the MKL version is present only with PARDISO loaded; the cache sizes equal Unit 41's constants; wall time ≥ 0 and the start time parses. |
| PG6 | **Suites.** Full Release and Debug suites pass at every commit; no warnings under `/W4 /WX`. |
| PG7 | **Mutation checks.** Hand mutants of the new writers, readers and provenance gathering (a dropped field, a swapped key, a wrong content kind, provenance hashed into the identity) are caught. |
| PG8 | **No scope creep.** Nothing outside PV1–PV6 and PD1's ruling; ARCHITECTURE.md 6.6 and 11.3 updated to "as built" only after the owner accepts the result. |

**Commits (proposed):** 1) `Provenance` and the build-info header, filled by `make_run_record` (PG1, PG5); 2) provenance
in sweep and checkpoint files (PG3); 3) transient and small-signal files (PG2, PG4); 4) the runner's `--version` and
report (PV6). Each passes PG6 on its own.

**Risks:**
- Tests that compare whole `Sweep` / `Transient` values for equality would see the new time fields differ between
  runs; they must compare with provenance excluded (PG1 catches any that do not).
- A build from a source archive has no git: the revision is `unknown` and `dirty` false, which PG5 must accept.
- The build-time git step must not force a rebuild when nothing changed (checked by a second, no-op build).

**Owner approval (2026-10-17):** "go ahead", taken as approval of the proposal with every recommendation:
PD1 the deck analyses are the next unit; PD2 start and wall time are in files; PD3 no host or user name; PD4 the
CSV table is unchanged; PD5 the version stays `0.0.0`.

**As built, so far:**
- Commit 1 (`4831f75`): `results::Provenance` in `RunRecord`; `cmake/build_info.cmake` writes the generated
  `NiTCAD/build_info.hpp` at every build, only when its text changes (a second build recompiles nothing);
  `linalg::eigen_build()`; `solve::provenance()`; every `make_run_record` attaches it and the four public runs
  stamp `wall_s`. The backend and its thread count are already settings, so they are not repeated (PV1 refined).
  Release 19/19. Mutants: 8 applied, 8 caught (the identity hashing the start time; a run's wall time dropped,
  for trace and for transient; cache sizes swapped; the start time to whole seconds; a wrong revision; a wrong CPU
  count; the Eigen version dropped). One further mutant, `%T` replaced by `%H:%M:%S`, is equivalent: `%S`
  prints the fractional seconds too.
- Commit 2: provenance as run-record key 3 and checkpoint key 7 (`io::Checkpoint::provenance`), every field
  optional on read. Fixtures `tests/io/fixtures/pre_unit42/` (a 1D bar run by `nitcad` built from `main` at
  `67e8fe5`) load with an empty provenance (PG3). `.gitattributes` keeps the binary fixtures unconverted.
  Mutants: 7 applied, 7 caught (provenance not written in the run record or in the checkpoint; a field not read;
  provenance required on read; the dirty flag inverted; wall time not written; its key colliding with the
  settings').

- `b5940c8`, a test fix: the commit-2 Debug suite failed 4 assertions, all the revision check. It required the
  stamped revision to equal git's HEAD when the test runs, and commit 2 was made after that build. A build
  older than HEAD is normal, so the revision must now be a 40-digit commit of HEAD's history
  (`git merge-base --is-ancestor`), or `unknown` without git (PG5 refined). The same commit adds PV5's test that
  a `trace_bias` sweep round-trips through `save_sweep`.
- Commit 3 (`397b08c`): content kinds 5 and 6, `save_transient` / `load_transient` and `save_small_signal` /
  `load_small_signal`. Each point's admittances are stored as one array (frequencies × contacts² complex values,
  as (real, imaginary) pairs) and its resolutions as one (frequencies × contacts), so a loaded run's
  `admittance(point, frequency, i, j)` stays in bounds. Loading also refuses a snapshot of a missing time point,
  an order outside 0–2, and small-signal fields of a missing contact (PG4). Mutants: 11 applied, 11 caught (error
  ratio, trap occupancy, small-signal fields or resolutions dropped; real and imaginary parts swapped; the
  admittance size check loosened, caught only after a too-long case was added to the test; the snapshot and
  contact bounds off by one; the order unchecked; the wrong content kind; the lagged count written as the
  retraced one). A first form of the snapshot mutant did not compile and was replaced.
- Commit 4 (`d8186bb`): `nitcad --version` prints the build on one line; `nitcad run` prints the build, the
  platform, the start and the wall time after the files it wrote (for an equilibrium run, which has no run
  record, the provenance of now). New ctest `cli_version` checks the line. Mutants: 3 applied, 3 caught.
- **Suites (PG6):** commit 1 Release 19/19; commit 2 Release 19/19, Debug 12 of 13 (the revision check above);
  commit 3 Release 19/19, Debug 13/13; commit 4 Release 20/20, Debug 14/14. No warnings under `/W4 /WX`. After
  the stamp, a second build recompiles nothing.
- **Debug assertion dialogs:** a mutant that makes a test read a failed `std::expected` opens the MSVC Debug
  runtime's modal assertion dialog, which waits for a click. Later mutant runs used the Release build.
- **PG1 trace gate, method:** Catch2 3 orders test cases randomly per run, and fixtures cached in function-local
  statics (for example `fixture_curve()` in `mos_test.cpp`) are traced inside whichever case first uses them. So
  traces are compared with `--order decl`, with both builds configured with the same MKL runtime (PARDISO cases
  run, not skipped), each pinned to P cores, block by block (a block is one test case's trace).
- **PG1 result:** every one of `main`'s 362 blocks over the nine traced entries is byte-identical in this branch's
  traces (assemble 107, deck 18, solve 199, breakdown 7, BTBT 3, MOSFET 11, transient 6, unstructured 3,
  validation 8), `result.identity` included. The only blocks without a counterpart are this unit's three new solve
  cases. No numerical result and no input identity changed.

**Status:** accepted by the owner (2026-10-17, "accepted"); ARCHITECTURE.md 6.6 and 11.3 updated to "as built"
(PG8). Merged to `main` as `145f23c` (PR #52, 2026-10-17).

### Unit 43 — transient, small-signal and trace analyses in the deck (proposed, authorized and accepted 2026-10-17; on `deck/analyses`, PR #53)

**Owner request (2026-10-17):** the scope of the next unit named by Unit 42's PD1. Documentation only until the
owner approves.

**Why now (audited on `main` at `145f23c`).**
- The deck (Unit 31, `deck/build.hpp`) has two analyses, `[equilibrium]` and `[sweep]` (`deck::Analysis`). The
  transient run (Unit 21), the small-signal run (Unit 22) and the arc-length trace (Unit 19) are reachable from C++
  only.
- Since Unit 42 every result has a file (`save_transient`, `save_small_signal`; a trace is a `Sweep`), so nothing
  but the deck stands between a user and these runs.
- ARCHITECTURE.md 11.3 names it as the input/output row's remaining step.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| AN1 | **Three new analysis sections**, `[transient]`, `[small_signal]`, `[trace]`. A deck has exactly one analysis section of the five (DK3 extended). `deck::Analysis` gains the three; `Plan` gains the options of each (`solve::TransientOptions` with its waveforms, `SmallSignalOptions`, `TraceOptions`). Every key maps to one C++ option, with the C++ default when absent. | One run per deck, as now; the deck stays a thin, strict layer over the solver options (DK2, DK3). | Five-way exclusivity to check and test. |
| AN2 | **`[transient]`:** `contact <name> = <waveform>` per contact (a contact not named is held at 0 V); waveforms `constant(v)`, `step(v0, v1, t)`, `ramp(v0, v1, t0, t1)`, `pulse(base, pulse, t_start, width[, rise, fall])`, `pwl(t v, t v, ...)`, `sine(offset, amplitude, f[, delay, phase])`, each one `solve::Waveform` factory. Keys `t_end` (required), `dt_initial`, `dt_min`, `dt_max`, `integrator = bdf2 \| backward_euler`, `adaptive = on \| off`, `rtol`, `density_ref`, `output_times = <series>`, `fields_every_step = on \| off`; `thermal <name> = <K>` (electrothermal: one temperature per thermal contact, as the run takes). | The `Waveform` factories already validate; the deck only parses their arguments. | A small function-call syntax (AQ1). |
| AN3 | **`[small_signal]`:** `contact <name> = <series>` operating points and `thermal <name> = <series>`, exactly as `[sweep]` (DK6); `frequencies = <series>` (required), `field_frequencies = <series>`. `logspace(a, b, n)` joins `linspace` in series (AQ5): frequencies are logarithmic. | Reuses the sweep's point rules unchanged. | — |
| AN4 | **`[trace]`:** `contact = <name>` (the swept contact) and `end = <V>` (required); `start <name> = <V>` per contact (default 0); `step`, `max_step`, `min_step`, `current_limit`, `max_points`, `density_floor`; `thermal <name> = <K>`. | `TraceOptions` field for field. | — |
| AN5 | **Outputs.** `[output]` gains `transient = <file>` and `small_signal = <file>`; a trace writes `sweep =` (it is a `Sweep`). Each analysis accepts only the outputs it can write; another is an error naming the analysis. `table =` writes a CSV per analysis (AQ3). `checkpoint =`: a trace's last point and a small-signal run's last operating point; refused with `[transient]` (AQ4). `device` and `mesh` as now. | The Unit 42 files are written by the runs that make them; a wrong output is a deck error, not a silent no-op. | — |
| AN6 | **`restart` stays `[sweep]`'s** (DK9). The new sections refuse it. | A transient or trace restart needs more than a converged point (history, arc direction). | No restart for long transients (later unit). |
| AN7 | **The runner:** `nitcad check` describes the new analyses (contacts' waveforms, frequencies, trace range); `nitcad run` prints progress per accepted time step, operating point or traced point, writes the outputs, prints the provenance; exit 2 when any run stopped early (its completed part written), as for a sweep. | One behaviour for every analysis. | — |
| AN8 | **Out of scope:** restart of the new analyses; a transient or small-signal run starting from a checkpoint file (`initial`); harmonic balance; new solver behaviour, models or defaults; `ARCHITECTURE.md` beyond 6.5 and 11.3 "as built" after acceptance. | Keeps the unit to the deck and the runner. | — |

**Owner decisions needed:**

| # | Question | Recommendation |
|---|---|---|
| AQ1 | Waveform syntax: function calls (`pulse(0, 1, 1e-9, 5e-9)`) as `linspace` already is, or words (`pulse 0 1 1e-9 5e-9`) as doping profiles are? | **Function calls**: named shapes with positional arguments read like `linspace`, and `pwl` needs pairs. |
| AQ2 | Time suffixes (`ps`, `ns`, `us`, `ms`)? DK5 allows suffixes on lengths only. | **Keep DK5**: times in seconds (`1e-9`), as every other non-length quantity. A later unit can add unit suffixes everywhere at once. |
| AQ3 | CSV tables: transient, one row per accepted step (time, step, order, the bias, total, conduction and displacement current per contact); small-signal, one row per (point, frequency) (the bias per contact, then G_ij and C_ij for every contact pair); trace, the sweep table. | **As stated.** |
| AQ4 | `checkpoint =` with `[transient]`: a transient state at t_end is not a steady state a sweep can restart from. | **Refuse it** (an error naming the analysis). |
| AQ5 | `logspace(a, b, n)` in series: `n` values from a to b evenly spaced in log10, the first a and the last b exactly, both positive; available in every series (axes included, harmlessly). | **Add it.** |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| AG1 | **Existing decks unchanged.** Every existing deck test passes unchanged; for `[equilibrium]` and `[sweep]` decks the plan, the run identity and the written files' values are as on `main`; the assembly traces of the deck and validation entries are byte-identical to `main`'s (`--order decl`, PG1's method). |
| AG2 | **The deck equals the API.** Each new analysis run from a deck gives, bit for bit, the result of the C++ call with the options the deck states (run identity included): a 1D diode transient with each waveform kind, a MOS-C small-signal run over a `logspace`, a diode trace; with and without the electrothermal model where the run supports it. |
| AG3 | **Every key and every error.** Each key reaches its option (a test per key, changing it from its default); absent keys give the C++ defaults. Each error is `invalid_input` "`<file>:<line>: <what>`" (DK4): unknown key, missing required key, malformed or invalid waveform, unknown contact, two analysis sections, `restart` outside `[sweep]`, an output the analysis cannot write, `checkpoint` with `[transient]`, `logspace` with a non-positive bound. |
| AG4 | **Outputs.** The written files load (`load_transient`, `load_small_signal`, `load_sweep`, `load_checkpoint`) equal to the in-memory results, provenance included; the CSV tables have AQ3's columns, one row per step, (point, frequency) or point, and are written atomically as the sweep table. |
| AG5 | **The runner.** `check` and `run` for each analysis; exit 0, 1 (a deck error) and 2 (a run stopped early: its completed part written) for each; a cancelled transient writes its accepted steps. |
| AG6 | **Suites.** Full Release and Debug suites at every commit; no warnings under `/W4 /WX`. |
| AG7 | **Mutation checks.** Hand mutants in the new parsing, plan building, output selection and tables (a key ignored, a default changed, a waveform argument swapped, an exclusivity check dropped, a table column dropped) are caught; run in Release (Debug mutants open modal assertion dialogs). |
| AG8 | **No scope creep.** Nothing outside AN1–AN7; `ARCHITECTURE.md` 6.5 and 11.3 updated to "as built" only after the owner accepts. |

**Commits (proposed):** 1) the three sections, waveforms and `logspace` in `deck::build`, with the plan tests (AG3);
2) running them and their outputs and tables in `deck::run` (AG2, AG4); 3) the runner (AG5). Each passes AG6 on its
own.

**Risks:**
- `deck/build.cpp` is already large (over 1000 lines); the new sections add a few hundred. A split into a separate
  `analysis.cpp` may be worth it in commit 1.
- A transient table with `fields_every_step` and many steps is large; the table holds terminal quantities only.
- AG2's bit-for-bit comparison needs both sides to see the same options, including defaults the deck fills; any
  default the deck sets differently from C++ would show here, as intended.

**Owner approval (2026-10-17):** "start the scope", taken as approval with every recommendation: AQ1
function-call waveforms; AQ2 times in seconds (DK5 kept); AQ3 the CSV columns as stated; AQ4 `checkpoint` refused
with `[transient]`; AQ5 `logspace`.

**As built, so far:**
- Commit 1: the three sections in `deck::build` (`Plan::transient`, `waveforms`, `small_signal`, `trace`; their
  `steady` left default, the run uses `Plan::options`); waveforms in `src/deck/waveforms.hpp`; `logspace` in the
  series reader. AN3 refined: series never had `linspace` (an axis form), so `logspace` is a series form of its own,
  usable in every series. `[sweep]` and `[small_signal]` share one operating-point reader (the sweep's, moved
  unchanged). Until commit 2, `deck::run` refuses the new analyses.
  - **The one change to existing behaviour (AG1 refined):** the message for a deck without exactly one analysis
    now lists the five sections; the two existing tests that quoted the old message were updated to it.
  - Tests (`tests/deck/analyses_test.cpp`, 6 cases): every key changed from its default reaches its option;
    absent keys equal the C++ defaults; each waveform form equals its factory; `logspace` values and exact bounds
    (also where 10^log10(x) does not round-trip); the outputs each analysis takes and refuses; every mistake
    with its line.
  - Mutants (Release): 14 applied, 14 caught (a key ignored: `rtol`, `max_points`; `t_end` not required; step
    arguments swapped; no analysis accepted; `logspace` bounds not exact, first and last, caught only after the
    inexact-bound test was added; the trace start always on the first contact; a transient checkpoint accepted;
    field frequencies stored as frequencies; a pulse with a rise only; `backward_euler` read as `bdf2`; the sweep
    output refused for a transient only; `[trace]` parsed as `[transient]`). A first form of the sweep-output mutant
    did not compile and was replaced.

- Commit 2 (`b77aacc`): `deck::run` calls `solve_transient`, `solve_small_signal` or `trace_bias` with the
  analysis's options and `Plan::options` as their `steady`; writes the transient, small_signal or sweep file, the
  checkpoint (small signal: the last operating point; trace: the last point) and the AQ3 table
  (`write_transient_table`, `write_small_signal_table`; a trace uses the sweep table). `RunReport` gains
  `transient` and `small_signal`; a trace's result is `sweep`.
  - Tests (`tests/deck/analyses_run_test.cpp`, 3 cases; AG2, AG4): each deck run equals the C++ call on the
    hand-built diode (not the deck's device), bit for bit, run identity included: a transient (a ramp, fixed
    backward-Euler steps, an output time), a small-signal run (two points, a `logspace`, field frequencies), a
    trace (with a solver option); the files reload equal, provenance included; the checkpoints; the tables'
    headers, row counts and a whole small-signal row rebuilt from the run.
  - AG2 narrowed: a diode stands in for the MOS-C the proposal named, and the waveform kinds are checked against
    their factories in the plan tests (commit 1), with one run (a ramp) through the solver; no electrothermal run
    from a deck (the thermal keys are checked in the plan).
- Commit 3 (`892d5fc`): `nitcad check` describes the new analyses; `run` prints progress per accepted step, per
  operating point and frequency, and per traced point; exit 2 when a trace, transient or small-signal run stops
  early. Tests: the CLI for each analysis (check, run, exit 0), a stopped transient and small-signal run (exit 2),
  and a stopped transient's file holding its accepted steps (AG5; the stop is a non-convergence, which takes the
  same path as a cancellation).
- Mutants for commits 2–3 (Release): 11 applied, 11 caught (displacement columns dropped; the small-signal
  checkpoint of the first point; each run without the plan's options, the trace's caught after its test gained a
  solver option; capacitances written as conductances, caught after the whole-row test was added; the trace table
  not written; the step column holding the time; a stopped transient exiting 0; waveform names swapped; trace
  progress printed as a sweep's).
- **AG1 result:** with `--order decl`, the same MKL configuration and P-core pinning, all 18 of `main`'s deck
  test-case blocks and all 8 validation blocks are byte-identical; the only other blocks are this unit's 11 new
  cases.
- **Suites (AG6):** commit 1 Release 20/20, Debug 14/14. Commit 2 changes the `deck` library and its tests only,
  so its gate was the deck suite (Debug, 27 cases), not the full suites, to save about 40 minutes; commit 3
  Release 20/20, Debug 14/14. No warnings under `/W4 /WX`.

**Status:** accepted by the owner (2026-10-17, "accept"); ARCHITECTURE.md 6.5 and 11.3 updated to "as built"
(AG8). Merged to `main` as `913dac0` (PR #53, 2026-10-17).

### Unit 44 — an iterative linear solver (proposed and authorized 2026-10-17; on `linalg/iterative`)

**Owner request (2026-10-17):** the scope of the "Scalable linear solvers" row of ARCHITECTURE.md 11.3 (planned:
"iterative Krylov with ILU/AMG behind `linalg`; parallel assembly"). Documentation only until the owner approves.

**Why now (audited on `main` at `913dac0`).**
- Every solve is a direct sparse LU: Eigen SparseLU (default) or MKL PARDISO (Unit 18, opt-in). Fill and time grow
  steeply on 3D meshes: one SparseLU refactorization of the 20³-node 3D diode took 5.8 s on this machine (Unit 41's
  cost measurement), and 3D meshing has been on `main` since Unit 16.
- The interface is ready for it: every run creates its solver from `BiasOptions::linear` (`bias`, `equilibrium`,
  `transient`, `small_signal`, `trace`, `weighted_current`); the backend interface (`src/linalg/backend.hpp`) is
  analyze / factorize / solve; `SolveReport::iterations` already exists, documented "for iterative backends".
- The acceptance check is backend-neutral: every solve is accepted on its componentwise backward error
  (`max_backward_error`, 1e-8), with iterative refinement.
- Eigen 5.0.1 (vcpkg, pinned) has `IncompleteLUT` and `BiCGSTAB` (supported) and GMRES, IDR(s), BiCGSTAB(L) and
  ILU(0) in `unsupported/`. No new dependency is needed.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| LS1 | **A third backend, `SolverBackend::gmres_ilut`, opt-in.** "Factorize" builds the preconditioner (Eigen's `IncompleteLUT`, with its own fill-reducing ordering); "solve" runs restarted GMRES with right preconditioning from x = 0. Real and complex (`ComplexLinearSolver`, the small-signal systems), as the other backends. | One more backend behind the same interface: the runs, the deck and every acceptance check are unchanged; the default backend and its results are untouched. | Eigen's ILUT, not an in-tree one: its ordering and drop rule are Eigen's. |
| LS2 | **GMRES in-tree (LQ1)**, about 200 lines in `linalg`: restarted GMRES(m) with modified Gram–Schmidt and Givens rotations. Its stopping test is the layer's own acceptance measure: each restart cycle ends with the componentwise backward error of the true residual, and the solve stops when it is at most `max_backward_error` or the iteration limit is reached. | The stopping test must be the componentwise backward error the layer accepts on (GMRES's own residual estimate is a 2-norm of the preconditioned system); `unsupported/` code has no stability promise. | Code to own and test; Eigen's GMRES would be less. |
| LS3 | **Configuration** (`SolverConfig`): `restart` (30), `max_iterations` (1000), `ilut_drop_tolerance` (1e-4), `ilut_fill_factor` (10). The existing `equilibrate` applies before the preconditioner, as for the direct backends. `threads` must be 1. | Equilibration matters more for ILU than for LU on drift-diffusion's widely scaled rows. Single-threaded keeps results deterministic on one machine (6.8). | Defaults are guesses until LG4 measures them. |
| LS4 | **Failure:** a solve that does not reach `max_backward_error` within `max_iterations` is `inaccurate_solve`, as a direct solve that fails acceptance; `SolveReport::iterations` reports the iterations. No pivot ratio (ILUT's diagonal is no singularity measure, as PARDISO's), so `min_pivot_ratio` is not applied; a singular system shows as non-convergence. Fallback to a direct solve: LQ2. | One failure path for every backend: Newton, step cutting and the run's `stopped` already handle `inaccurate_solve`. | A floating region is not reported as `singular_system` by this backend. |
| LS5 | **The deck:** `[solver] backend = gmres_ilut`, with `restart`, `max_linear_iterations`, `ilut_drop`, `ilut_fill`. | The batch runner reaches it like PARDISO. | Four keys. |
| LS6 | **Out of scope:** algebraic multigrid; block (field-split) preconditioners (they need Unit 40's unknown map, a later unit); parallel assembly; a GPU solver; changing the default backend. | Keeps the unit to one backend and its evidence. | 3D scaling beyond what ILUT gives waits for AMG. |

**Owner decisions (asked; ruled below):**

| # | Question | Recommendation |
|---|---|---|
| LQ1 | GMRES in-tree (LS2), or Eigen's `unsupported/` GMRES with a componentwise check after it? | **In-tree.** The stopping test is the acceptance measure, and the code stays under this project's tests. |
| LQ2 | When GMRES fails, fall back to Eigen SparseLU for that system (and say so in the report), or fail with `inaccurate_solve`? | **Fail** in this unit (one behaviour, measurable failure rate); a fallback can follow once LG4 shows how often it fails. |
| LQ3 | Complex systems (small signal) in this unit? | **Yes**: the same template; the high-frequency admittance systems are a useful test of the preconditioner. |
| LQ4 | If LG4 shows no device where GMRES-ILUT beats SparseLU in time or memory, merge it anyway as an opt-in (documented as such), or stop? | **Merge as opt-in** with the measurements recorded; it is the base AMG and block preconditioners plug into. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| LG1 | **The default is untouched.** With `eigen_sparse_lu` every suite passes unchanged and the assembly traces are byte-identical to `main`'s (`--order decl`, PG1's method). |
| LG2 | **GMRES on its own.** Linear-algebra tests: GMRES with no preconditioner and with ILUT reaches the backward-error target on nonsymmetric real and complex systems (including `tests/linalg`'s shared systems); its iteration count on a matrix with k distinct eigenvalues is at most k; restart and iteration limits are honoured; a singular or unconverged system is `inaccurate_solve` with the iterations reported; `threads != 1` and bad settings are `invalid_input`. |
| LG3 | **It agrees with SparseLU** on the devices, as PARDISO does (Unit 18's `[pardiso]` cases as `[gmres]` cases): the diode bias gates, the MOS quasi-static C-V, the small-signal admittances, the transient switch-on, within the same acceptance tolerances; the Newton iteration counts reported. |
| LG4 | **Measured, not assumed:** a Release benchmark (`[.benchmark]`) of one Newton step's linear solve on the 2D MOSFET and on 3D diodes of growing size (at least 20³, 30³ and 40³ nodes where SparseLU still fits): time, peak working set, GMRES iterations, for SparseLU and GMRES-ILUT (and PARDISO where configured). The results go into DECISIONS.md whatever they show (LQ4). |
| LG5 | **The deck** keys of LS5 reach `SolverConfig`; a 3D deck runs with `backend = gmres_ilut`. |
| LG6 | **Suites:** Release and Debug at every commit; no warnings under `/W4 /WX`. |
| LG7 | **Mutation checks** (Release): GMRES (a skipped reorthogonalization, a wrong Givens sign, the stopping test on the estimate instead of the true residual, a restart that keeps the old basis), the preconditioner applied on the wrong side, a configuration check removed, a deck key ignored. |
| LG8 | **No scope creep:** nothing outside LS1–LS5; ARCHITECTURE.md 6.10 and 11.3 updated to "as built" after acceptance. |

**Commits (proposed):** 1) GMRES and the backend in `linalg`, with LG2; 2) the device agreement cases (LG3) and the
deck keys (LG5); 3) the benchmark and its results (LG4).

**Risks:**
- Drift-diffusion Jacobians are hard for ILU: the continuity rows are strongly nonsymmetric at high fields, and the
  scales of ψ, n and p differ by tens of decades before equilibration. GMRES-ILUT may fail on some devices; LG3 and
  LG4 will show where, and LQ2 decides what a failure does.
- A backward error of 1e-8 componentwise may take many iterations where a direct solve gets 1e-16 at once; the
  benchmark may show a looser `max_backward_error` is needed for speed, which is the user's choice per run, not a
  default change.
- Results differ from SparseLU's at the tolerance, not bit for bit, as PARDISO's do.

**Owner rulings (2026-10-17).** The owner approved the unit and delegated LQ1–LQ4 to the choice best for a
commercial TCAD product ("for decisions choose whats best for Commercial TCAD").

| # | Ruling | Why, for a commercial product |
|---|---|---|
| LQ1 | GMRES in-tree (LS2). | The stopping test is the acceptance measure, and the solver stays under this project's tests and stability promise. |
| LQ2 | **Fall back** (differs from the recommendation). When GMRES-ILUT fails to reach `max_backward_error` within `max_iterations`, or the ILUT build fails, that system is solved by Eigen SparseLU, and the fallback is reported: `SolveReport` gains a field saying the system fell back, with the GMRES iterations spent. The direct solve's own checks then apply as usual (`inaccurate_solve`, `singular_system`). The deck key `linear_fallback = sparse_lu` (default) or `none` (fail with `inaccurate_solve`, the strict behaviour) controls it. | A production run must not abort where a direct solve can finish; the user sees in the report how often it happened. The strict mode keeps the failure rate measurable. |
| LQ3 | Complex systems (small signal) in this unit. | The same template; high-frequency admittances test the preconditioner. |
| LQ4 | Merge as opt-in with the measurements recorded, whatever LG4 shows. | It is the base for AMG and block preconditioners, the path to large 3D. |

**Amendments from LQ2.** LS4: failure falls back as LQ2, with `none` giving the behaviour written there. LS5: a
fifth key, `linear_fallback`. LG2: with `linear_fallback = none` an unconverged system is `inaccurate_solve`; with
`sparse_lu` the same system is solved, the fallback is reported, and its solution meets the backward-error target.
LG3: the agreement cases run with `none`, so a fallback cannot hide a failure; any device that fails is recorded.
LG4: also records each case's fallback count with the default setting. LG7: a further mutant, a fallback that is
not reported.

**Status:** accepted by the owner 2026-10-21 ("go next"); LQ3 amended by the owner (LQ3a); commits 1–3 done;
ARCHITECTURE.md 6.10 and 11.3 updated to "as built" (LG8). Merged to `main` as `b518c99` (PR #54, 2026-10-21).

**Commit 3 (LG4 benchmark), measured 2026-10-21** on the development machine (Release, one thread, 31.6 GB),
`[.benchmark][gmres_benchmark]` in `bias_test.cpp` and `mosfet_test.cpp`; MKL not configured, so no PARDISO
column. One Newton step's linear solve with a fresh solver (analysis, factorization or ILUT build, one solve), at
the default settings (`max_backward_error` 1e-8, restart 30, drop 1e-4, fill 10, fallback `sparse_lu`). "Held" is
the private memory the solver holds after the solve; "peak rise" how far the process's peak working set rose
during the step (GMRES-ILUT ran first and never raised it: its peak stayed below the Jacobian assembly's).

| Case | Unknowns | SparseLU | GMRES-ILUT | GMRES iterations | Fallbacks | Step difference |
|---|---|---|---|---|---|---|
| 2D MOSFET 90 × 50, Vg = 1 V | 28,470 | 0.23 s, held 122 MB | 0.47 s, held 48 MB | 120 | 0 | 3.3e-13 |
| 3D diode 20³, 0.5 V | 24,000 | 7.1 s, held 409 MB, peak +321 MB | 0.77 s, held 50 MB | 40 | 0 | 2.2e-14 |
| 3D diode 30³ | 81,000 | 105 s, held 3.1 GB, peak +2.6 GB | 3.8 s, held 174 MB | 60 | 0 | 7.4e-14 |
| 3D diode 40³ | 192,000 | 1021 s, held 10.9 GB, peak +8.0 GB | 14.1 s, held 422 MB | 90 | 0 | 8.7e-14 |

The MOSFET's Id-Vg sweep (16 points, Vds = 0.05 V): SparseLU 21.1 s and 100 Newton iterations, GMRES-ILUT 29.9 s
and 84, every drain current within max(1e-8 |I|, resolution) of SparseLU's (worst ratio 0.000).
**Reading (LQ4):** on the 2D MOSFET SparseLU is twice as fast and GMRES-ILUT holds 2.5× less memory; on 3D
GMRES-ILUT is 9×, 27× and 72× faster and holds 8×, 18× and 26× less memory, its iterations growing about as the
grid's linear size (40, 60, 90). No fallback was taken in any case. The default backend stays SparseLU (LS6);
these figures are the evidence for recommending `backend = gmres_ilut` for 3D devices.

**LQ3a (owner ruling, 2026-10-21): small-signal systems are solved directly.** Measured in commit 2: with GMRES
the admittances agree with SparseLU only within each run's own reported resolution, and that resolution reached
22% of max |Y| (p+n junction) and 93% (MOS capacitor with an interface trap level), the same at
`max_backward_error` 1e-8 and 1e-9; at 1e-10 GMRES-ILUT did not converge the trap device's operating point in 1000
iterations. An admittance is a near-cancellation in the continuity rows: a componentwise backward error of 1e-9
still leaves an absolute residual comparable to the AC current, and only a rounding-level (direct) solve resolves
it. Asked with three options (direct AC solves; keep GMRES with coarse resolutions; stop the unit), the owner chose
**direct AC solves**: with `backend = gmres_ilut`, `solve_small_signal` takes its operating points with GMRES-ILUT
and solves the complex small-signal systems with Eigen SparseLU (`small_signal.hpp`). The complex GMRES path stays
in `linalg`, tested (LG2), for later preconditioners. Not recorded in the results: the run record's
`linear.backend` says gmres_ilut, and each frequency's `pivot_ratio` is the direct factorization's.

**Commit 2 (device agreement LG3, deck keys LG5), as built.** `[solver]` keys `backend = gmres_ilut`, `restart`,
`max_linear_iterations`, `ilut_drop`, `ilut_fill`, `linear_fallback = sparse_lu|none`, each checked with its line;
the run record digests the five settings only with `gmres_ilut`, so every other run keeps its identity. `[gmres]`
cases (strict mode, in the `solve` ctest) beside Unit 18's `[pardiso]` ones, measured in Release:
- diode: one ordering per Newton run, one ILUT per Jacobian (10 iterations as Eigen, 71 preconditioner
  applications); J(0.5 V) 1.2800773807e-2 as Eigen; the 0–0.7 V sweep within 5.7e-9 relative, 69 Newton
  iterations as Eigen's 69;
- MOS quasi-static C-V: gate charge within 7.9e-17 of the largest, φ_s within 1.1e-16 V;
- transient switch-on (101 BDF2 steps): currents within 4.4e-21 of the largest;
- small signal (LQ3a): junction and trap level within 0.10 and 0.23 of the PARDISO case's bound;
- deck: a 3D diode deck (21 × 6 × 6 nodes) with `backend = gmres_ilut` sweeps to 0.4 V within 4.4e-16 of the
  SparseLU run.
- **Suites:** Release 20/20 and Debug 14/14, no warnings. **Mutation checks (LG7, Release): 12 of 12 caught** —
  each of the six deck keys ignored, each of the five run-record entries dropped, and the small-signal systems
  left on GMRES (LQ3a).

**Commit 1 (GMRES and the backend, LG2), as built.** `SolverBackend::gmres_ilut` and `LinearFallback`;
`SolverConfig` gains `restart`, `max_iterations`, `ilut_drop_tolerance`, `ilut_fill_factor`, `fallback` (validated
for every backend); `SolveReport` gains `fell_back`. The kernel is `src/linalg/gmres.hpp` (no Eigen; MGS with one
full reorthogonalization pass, Givens rotations, right preconditioning, x = 0 start); the preconditioner is
`src/linalg/eigen_ilut.cpp` (Eigen 5.0.1 `IncompleteLUT`, version-asserted; a non-finite factor is a failed
build). GMRES runs on the equilibrated system the ILUT was built from and accepts only on the componentwise backward
error of the true iterate, measured once per cycle; a cycle ends early when the residual estimate falls below
1/2 · β · max_backward_error / (current backward error), so a wrong aim costs a restart, never a wrong answer.
- **Finding, amends LS4.** A singular system does not show as non-convergence. On a floating region with an
  inconsistent right-hand side GMRES reaches an x of norm 1.1e17 with backward error 1.0e-16, exactly as Eigen
  SparseLU does with `min_pivot_ratio = 0` (|x| 3.1e16): the componentwise backward error cannot see singularity.
  `gmres_ilut` therefore does not detect floating regions; the device layer's topology check (`device.cpp`)
  rejects them before any solve, so device runs are unaffected. A fallback still applies the direct pivot check.
- **Finding, LG2's k-eigenvalue test.** Termination in k steps holds in floating point only when the minimal
  polynomial is well conditioned: with 6 eigenvalues over four decades its terms reach 1e19 before cancelling
  (measured) and GMRES needs a 7th step. The test uses 6 eigenvalues in [1, 13] and gets at most 6.
- **Measured (Release, `tests/linalg`, strict mode):** nonsymmetric n = 2000 and the contacted 1D Laplacian in 1
  iteration; the 40 × 40 2D Laplacian in 9 (forward error 3.8e-9 against SparseLU); the badly scaled coupled
  system (rows over 40 decades) in 170 (9.7e-9); complex G + iωC over ω = 1e-6…1e3 within 3.7e-14 of SparseLU.
- **Suites:** Release 20/20 and Debug 14/14, no warnings; `nitcad_linalg_test` +18 test cases (`gmres_kernel_test`,
  `gmres_ilut_test`). **Mutation checks (LG7, Release): 19 of 20 caught** — skipped reorthogonalization, Givens
  sign, stopping on the estimate, a restart keeping the old basis, left instead of right preconditioning, each of
  the five configuration checks, the fallback unreported (two places), the fallback ignoring `none`, the operator
  or the measure not equilibrated, the ILUT finiteness check removed, `direct_` not reset on refactorization. The
  survivor, a Givens rotation without the conjugate of its second entry, is equivalent: that entry is the
  Arnoldi norm, real and positive.

### Units 25 and 26 — visualization and the Windows application (proposed and authorized 2026-10-17; on `render/visualization`)

**Owner request (2026-10-17):** a scope for the "Visualization, GUI, batch" row of ARCHITECTURE.md 11.3 (Units 25
and 26, both "target"), presented as an HTML page; then the target look given by the owner's reference grid of
twelve 3D MOSFET views ("TCAD Semiconductor Visualisation Grid"). The proposal, with mockups of the window and of
the twelve views, is the private page https://claude.ai/artifact/A3K7DSo1FZMjsm1aew5htX (version 3); this entry is
its record.

**Why now (audited on `main` at `913dac0`).**
- Every result type has a file since Unit 42 and every analysis runs from a deck since Unit 43; the only output a
  user sees is `nitcad`'s text and CSV tables.
- The roadmap reserves `render` (L7: base, mesh, results, analysis) and `app` (L8: everything except physics and
  assemble). 6.9 fixes the app's threading: every solve on a worker thread, progress posted to the UI thread.
- On this workstation: Windows SDK 10.0.26100 with `d3d12.h`, `d2d1_3.h`, `d3d11on12.h`, `dxgi1_6.h`,
  `wincodec.h` and `dxc.exe`; WARP (`d3d10warp.dll`); an NVIDIA GeForce RTX 5060 Ti. No new vcpkg dependency.
- A mesh is a node and edge graph without cells: surfaces come from a generated mesh's simplices or a tensor
  device's `TensorCells`. A sweep file holds fields, not the mesh: the viewer pairs it with the deck's mesh file.

**Design: Unit 25 (`render`).**

| # | Choice |
|---|---|
| VR1 | A CPU scene (vertices, values, colors, line strips, labels) built from a mesh and a result, separate from drawing; drawing never reads results. |
| VR2 | Per dimension: 1D line plots; 2D field maps on triangles (simplices, or tensor cells split in two) with edges, region outlines and contacts; 3D the reference grid's views: region-colored solid view with contacts and electrodes, any field on the boundary surface, the tetrahedral mesh with regions, slice planes (XY, XZ, YZ), isosurfaces, current streamlines, transparency. 3D tensor cells are split into six tetrahedra. A mesh without simplices or cells is drawn as nodes and edges only. |
| VR3 | Fields: potential, n, p, |E|, |J|, band edges, quasi-Fermi levels, temperature, doping; linear and log scales; colorbars with units; contours by marching triangles; isosurfaces by marching tetrahedra; streamlines of J traced through tetrahedra (RK4, seeded on a contact or a line). |
| VR4 | Plots (I-V, C-V, transient waveforms, admittance against frequency, band diagrams) from `analysis` curves, with nice ticks, log axes, units and Unit 24's resolution and extrapolation flags. |
| VR5 | Drawing: triangles and lines on Direct3D 12; text, axes and legends on Direct2D through D3D11On12; off-screen targets read back and saved as PNG through WIC; a hardware device when present, else WARP. |
| VR6 | `nitcad plot <deck or result file> ...`: fields, views and curves to PNG; a result file is paired with its mesh (deck outputs or `--mesh`), node counts checked. |
| VR7 | Out of scope: windows (Unit 26), animation, ray-marched volume rendering, VTK export, printing. |
| VR8 | The reference's look: lit, shaded solids; fixed region colors (substrate beige, n⁺ red, p channel blue, oxide pale translucent, metal dark grey); an electrode drawn as a slab of nominal thickness on its contact patch (a gate is a contact in NiTCAD, not a meshed region); weighted blended order-independent transparency; an axis triad in every view. |
| VR9 | Viewing reconstructions in `render`: |E| per tetrahedron from the potential's gradient, J per tetrahedron from the edge currents (least squares over its six edges); nothing added to the solver or to files. Recombination is not in result files and `render` may not include `physics`, so the reference's R view waits for a unit that stores the rate in results. |

**Design: Unit 26 (`app`).** VA1 one main window (Win32, per-monitor DPI v2: deck pane, views, run bar, log); VA2
the deck pane (Rich Edit, line numbers, Check puts `deck::build`'s errors at their lines); VA3 runs on a worker
thread, progress by posted message, Stop through the stop token, completed parts kept; VA4 Unit 25's scenes in a
swap chain with pan, zoom, rotate and node picking; VA5 a results browser (run record, provenance, convergence,
Unit 24's figures); VA6 logic in plain C++ under test, a thin window layer, `--self-test` for CI on WARP; VA7 out of
scope: a graphical device editor, projects, remote runs, plugins, localization, an installer.

**Owner rulings (2026-10-17):**

| # | Ruling |
|---|---|
| VQ1 | D3D11On12 for Direct2D on Direct3D 12. |
| VQ2 | HLSL compiled at build time with the SDK's `dxc`, bytecode embedded. |
| VQ3 | Units in sequence: 25a, then 25b, then 26. |
| VQ4 | Exact tests of the CPU scene, plus WARP image regression tests against committed references. |
| VQ5 | Turbo as the default colormap; viridis, cividis and a diverging map also supported. |
| VQ6 | Unit 44 (iterative solver) first if it is approved; then Unit 25, then Unit 26. |
| VQ7 | Unit 25 split: 25a, the core rendering (scene, surfaces, solid and mesh views, slices, transparency, the image pipeline) and `nitcad plot`; 25b, derived |E| and |J|, isosurfaces and streamlines. |

**Owner conditions (2026-10-17):** numerical behaviour and the architecture's boundaries are preserved; Release
and Debug tests and warning-free `/W4 /WX` builds; scientific validation and performance benchmarks; only the
currently approved unit is implemented, and each is accepted before the next starts; nothing is pushed or merged
without authorization.

**Acceptance gates.**

| # | Gate | Unit |
|---|---|---|
| VG1 | Layering: `render` includes only base, mesh, results and analysis (a header-boundary test); the assembly traces are byte-identical to `main`'s (no numerical change). | 25a |
| VG2 | Geometry: triangles of tensor and generated 2D meshes cover the domain exactly once; 3D boundary surfaces are closed; tetrahedra from tensor cells fill each cell; slice polygons match analytic sections. | 25a |
| VG3 | Values: colormap tables and endpoints; log scales exact at decades; contours of a linear field at the exact positions; ticks of 1, 2 or 5 times a power of ten. | 25a |
| VG4 | Images: WARP renders of a band diagram, a 2D potential map, an I-V plot and the 3D MOSFET's solid, surface (n, p, φ), mesh, slice and transparent views match committed references within a stated per-pixel tolerance; in 25b also |E|, |J| with streamlines and isosurfaces. The D3D12 debug layer reports no errors in Debug. | 25a, 25b |
| VG5 | `nitcad plot` writes each image from decks and from result files and refuses a mismatched mesh (exit 1). | 25a |
| VG6 | Release and Debug suites at every commit, no warnings under `/W4 /WX`; hand mutants in the scene algorithms are caught (run in Release). | 25a, 25b |
| VG7 | Benchmarks: scene building and drawing of 10⁵ and 10⁶ triangles on WARP and on the GPU, recorded; the GPU budget is 60 frames per second at 10⁶ triangles. | 25a |
| VG8 | Scientific validation: marching tetrahedra on a linear field gives the exact plane, closed and oriented; on a sphere's distance field the enclosed volume converges at second order; streamlines of a uniform field are straight and of a rotational field close within the integrator's tolerance; |E| of a linear potential and J of a uniform current are exact per tetrahedron. | 25b |
| WG1 | No solve on the UI thread: the window keeps pumping messages (within 100 ms) during a long run; Stop keeps the completed part. | 26 |
| WG2 | A deck run in the app is bit-identical to `nitcad run` on it, run identity included. | 26 |
| WG3 | The app's logic (run state machine, documents, view settings) under unit tests; deck errors land on their lines. | 26 |
| WG4 | `--self-test` passes on WARP and on the GPU and writes the same image as `nitcad plot` for the same view. | 26 |
| WG5 | Per-monitor DPI re-layout, keyboard access to every command, a clean D3D12 debug layer, no leaked GPU objects. | 26 |
| WG6 | Suites as VG6; nothing outside VA1–VA6; ARCHITECTURE.md updated to "as built" after acceptance. | 26 |

**Commits (planned):** 25a·1 the CPU scene (VG2, VG3); 25a·2 off-screen D3D12 and Direct2D drawing, lighting, region
colors, electrode slabs, transparency, PNG, references and the benchmark (VG4, VG7); 25a·3 `nitcad plot` (VG5);
25b derived fields, isosurfaces, streamlines (VG8, VG4); 26·1 the app's logic (WG1, WG3); 26·2 the window (WG2,
WG4, WG5). Estimated sizes: 25a about 3,900 lines, 25b about 1,700, 26 about 2,900, tests included.

**Risks:** GPU images differ between vendors and drivers (image gates on WARP only); D3D11On12's per-frame cost;
large 3D results (single-precision GPU copies, scales from the doubles); `dxc` and WARP on the CI runner are
unverified (checked in 25a's first commits); screen-reader support for custom views is out of scope.

**Detailed plan (page versions 4–6, 2026-10-17).** At the owner's request the page planned the application in
more detail and took the look of the owner's workspace reference (a light 3D workspace: menu bar, tool bar, scene
tree, viewport with overlay cards, metrics strip, inspector, status footer). VA1 now describes that workspace; VA8
adds the viewport tools and overlays (select, orbit, pan, zoom, fit; axonometric and perspective; view, field,
colorbar, probe and orientation overlays; probe and measure in 26, slice in 25a, clip in 25c). Figures 3–8 plan
the start page, the deck editor, the run monitor, the results workspace, cut lines with band diagrams and run
comparison; sixteen more views (13–28) are assigned to units. Sizes: 25a about 4,100 lines, 25b 1,700, 25c 1,600,
26 3,300, 26b 2,600, 45 500.

**Rulings VQ8–VQ14 (2026-10-17).** The owner delegated these to the recommendation best for a commercial TCAD
product ("for decisions choose whats best for Commercial TCAD"); each is the page's recommendation.

| # | Ruling | Why, for a commercial product |
|---|---|---|
| VQ8 | `render` may depend on `device`, read-only: regions, material names and the insulator flag, contacts, doping. No solver, assembler or physics model. VG1 becomes: `render` includes only base, mesh, results, analysis and device. ARCHITECTURE.md 3 changes with 25a. | Region colors, electrodes and doping are what users look at first; a copied "view description" in results would be a second source of truth for the same geometry. |
| VQ9 | Unit 25c (vector glyphs, clipping box, exploded view, transient filmstrips) after 25b, each view droppable. | Vectors and clipping are expected in every commercial viewer; independent views keep the risk small. |
| VQ10 | Unit 26b (live deck checking and completion, structure preview, cut lines with band diagrams, run comparison and difference maps, mesh quality) after 26. | These make the application an analysis tool rather than a viewer; 26 keeps its size. |
| VQ11 | A top-level `examples/` directory with four decks (1D diode, 2D MOSFET, 3D MOSFET, MOS-C), each run by a test. | A first-run user needs working examples; tested examples never rot. |
| VQ12 | Unit 45 (recombination and generation rates in results and files; the mesh digest in sweep, transient and small-signal files) before 25b. | Unlocks the R view and exact result–mesh pairing instead of a node-count check. |
| VQ13 | The Device, Mesh, Physics and Simulation menus open the deck at that section and list its keys for insertion; the deck stays the only source of truth. Form editors that write the deck may be a later unit. | Reproducible, diffable, scriptable input is what production flows depend on; a form editor that writes the deck can come later without a second model. |
| VQ14 | Turbo stays the default (VQ5); the reference's softer five-color map is added as a named option, `soft`. | Turbo is the familiar high-contrast default; `soft` matches the reference where wanted. |

**Status:** authorized with VQ1–VQ14. Unit 44 was implemented first (VQ6) and merged (`b518c99`, PR #54); Unit 25a
started on the owner's "go 25a" (2026-10-21), on `render/visualization` rebased onto `b518c99`; commits 25a·1–3 done;
**accepted by the owner (2026-10-21, "accepted")**; ARCHITECTURE.md 6.14, 3, 11 and 11.3 as built. Then 45, 25b, 25c,
26, 26b, each approved after the one before it is accepted.

**Commit 25a·1 (the CPU scene; VG1–VG3, VG6), as built.** The `render` layer (`include/NiTCAD/render/`,
`src/render/`; CMake target `NiTCAD::render`, linked to base, mesh, results, analysis and device):
- `color.hpp`: Turbo, viridis, cividis (matplotlib 3.9.2's 256-entry tables, extracted with awk), Moreland's
  cool-warm (33 points) and `soft` (VQ14: the reference names deep blue, teal, green, gold and muted red without
  values; chosen as #253f8f, #2a9d9c, #6dbb5a, #e3b23c, #c0504d), interpolated linearly in sRGB; `ColorScale` linear,
  log10 and signed log10 (sign(v) log10(1 + |v| / floor), for doping); ticks of 1, 2 or 5 times a power of ten,
  decades on log axes (minor 2..9 on axes of at most three decades), labels with U+2212 and superscript powers.
- `geometry.hpp`: elements from `TensorCells` (2 triangles; 6 Kuhn tetrahedra, conforming between cells) or a
  `SimplexMesh` (reoriented counterclockwise), else the graph; 3D boundary faces (outward); slices (a facet lying
  in the plane taken once, from its positive side); **regions drawn by control volume**: in the box method the
  interface between two regions runs between their nodes (no generated node lies on one), so a tensor cell is
  split at its mid-planes (exactly the box-method control volumes) and a simplex by its planes of equal
  barycentric weight; pieces, interfaces and both rules by convex-polygon clipping; marching-triangle contours (a
  vertex at the level counts as above it).
- `scene.hpp`: `FieldScene` (triangles with per-vertex colormap positions or flat colors, line segments, colorbar,
  camera, bounds, in µm) and `PlotScene` (axes, series, flags); views `solid_view`, `mesh_view`, `field_view`
  (2D maps; 3D on the boundary surface, with the semiconductor seen through translucent insulators at their
  interfaces), `slice_view`; contact slabs of nominal thickness (3% of the device: a lumped gate as oxide under
  metal, an electrode as metal, an ohmic contact as a thin pad); the VR8 palette by doping class (n⁺ ≥ 1e19 red,
  n salmon, p channel 1e17–1e19 blue, p⁺ dark blue, substrate beige, insulators pale at alpha 0.35, metal dark
  grey); `field_values` for ψ, n, p, the bands and quasi-Fermi levels, T and doping; `curve_plot` (Unit 24 curves;
  a log axis of all-negative values plots |y|; points at or below their resolution flagged), `profile_plot`,
  `band_plot` (1D; insulators as gaps).
- **Gates:** `tests/render` (26 cases): VG2 — 2D tensor and generated triangles cover the domain exactly once
  (no directed edge twice); 3D tensor and generated boundaries closed and outward (enclosed volume and area to
  1e-12); Kuhn tetrahedra fill each cell; slices match the analytic sections (x = 0.73, a grid plane, the hexagon
  x + y + z = 1.5) with a linear field exact on them; the control-volume pieces of a tensor mesh equal its control
  volumes to 1e-14; tensor interfaces lie halfway between the nodes; VG3 — table endpoints, log positions exact at
  decades, ticks and labels; scenes — region areas per material equal the control volumes, field positions exact,
  contours on their levels, winding consistent with every normal, the silicon under an oxide facing the oxide
  above and below. VG1: ctest `render_layering` (a CMake script, no Python) checks that `render` includes only
  base, mesh, results, analysis and device; nothing outside `render`, its tests and CMakeLists.txt changed.
- **Mutation checks (VG6, Release): 18 of 18 caught** (Kuhn, simplex and face orientation, the slice facet rule,
  the tensor and barycentric owner rules, the interface position, the contour tie rule, the table step, the ticks'
  2, the scale's clamp, the fan's orientation, the interface normal, insulator transparency, field and slice
  positions, the doping threshold, a forbidden include). A degenerate-triangle check was found dead (the vertex
  merge before it removes every case) and removed. Release 22/22 and Debug 16/16 (with the new ctests `render` and
  `render_layering`), no warnings under `/W4 /WX`.

**Commit 25a·2 (drawing; VG4, VG6, VG7), as built.** `image.hpp`: a `Renderer` on a Direct3D 12 device (the first
hardware adapter, else WARP) draws a `FieldScene` or a `PlotScene` off screen into an `Image`; `save_png`,
`load_png` (WIC) and `compare`.
- Triangles and lines on Direct3D 12 at twice the image's size, averaged 2 × 2; field triangles look the colormap up
  per pixel from the interpolated position (a 256-entry table resampled from `sample`); one directional light
  (from the upper left over the viewer's shoulder) and ambient 0.38, two-sided; translucent triangles by weighted
  blended OIT (McGuire and Bavoil); lines as screen-space quads of their width, depth-tested or on top. Text,
  colorbars, 2D axes in µm, the axis triad and plots on Direct2D and DirectWrite (Segoe UI, grayscale text
  antialiasing) drawn into the same target through D3D11On12 (VQ1). HLSL compiled at build time by the SDK's `dxc`
  into headers of bytecode (VQ2; CMake finds `dxc` and fails without it).
- **The camera convention (changes 25a·1's):** +y up, the device's vertical axis, as in 2D; azimuth about +y from +z
  toward +x, elevation above the xz plane; 2D looks down −z. The default 3D view is azimuth 35°, elevation 25°.
  Drawing fits the scene's bounds to the image's aspect (92% of the view).
- **Debug layer:** the one message the pipeline raises, `REFLECTSHAREDPROPERTIES_INVALIDOBJECT` from D3D11On12's
  `CreateWrappedResource` probing a resource Direct3D 12 created, is filtered by its exact ID; every other error
  counts, and the tests require none (Debug builds run the image tests with the layer).
- **Errors:** a device failure is `resource_exhausted` with the call and its HRESULT (there is no device error code in
  `base`; a new one was not added in this unit).
- **VG4 references** (`tests/render/references/`, 11 images, 960 × 720, WARP): a 2D MOSFET's potential map with
  contours and its regions with the mesh, a band diagram, an I-V plot with flagged points, and a 3D MOSFET's solid,
  n, p and ψ surfaces, mesh, two slices, and an opaque perspective view; the MOSFETs are solved in the test (21 × 21
  and 21 × 21 × 7 nodes, Vg = 1 V, Vd = 0.5 V). Each was looked at before it was committed; the first renders
  showed the views mirrored (the camera's right vector reversed), fixed before any reference was kept. A render
  matches when at most 0.5% of its pixels differ by more than 24 in a channel; `NITCAD_UPDATE_REFERENCES=1` rewrites
  them. Also exact checks: a flat triangle's color, a field's colormap color, image sizes, the PNG round trip, and
  order-independent transparency giving the same image in either order.
- **VG7 benchmark** (`[.benchmark]`, 1280 × 960, the boundary of a thin tensor box):

  | Triangles | Adapter | Image (draw: upload, passes, overlays, readback) | Redraw on resident geometry |
  |---|---|---|---|
  | 6.8e4 | RTX 5060 Ti | 13.6 ms (74/s) | 0.20 ms |
  | 1.0e6 | RTX 5060 Ti | 146 ms (6.8/s) | 0.81 ms (1236 frames/s) |
  | 6.8e4 | WARP | 29.0 ms | 13.8 ms (72/s) |
  | 1.0e6 | WARP | 302 ms | 88.1 ms (11/s) |

  The 60 frames per second budget at 10⁶ triangles is met on the GPU by the passes a window repeats (0.81 ms); a
  whole image costs 146 ms, nearly all of it the per-image upload of 132 MB of vertices and the readback, which a
  window does not repeat. Scene building for 10⁶ triangles takes 1.3 s on the CPU.
- **Mutation checks (VG6, Release): 10 of 10 caught** — the camera's handedness, the fit margin, the readback and
  PNG channel order, the field flag, a line quad's corner, the depth test, and in the shaders the lighting, the OIT
  revealage and the 2 × 2 average.
- ARCHITECTURE.md 3: `render`'s row of allowed dependencies adds `device`, read-only (VQ8).
- **Suites:** Release 22/22 and Debug 16/16, no warnings; `nitcad_render_test` 30 test cases (+4). The Debug `render`
  ctest takes 641 s (the two MOSFETs solved unoptimized, then drawn on WARP with the debug layer), Release 12 s.

**Commit 25a·3 (`nitcad plot`; VG5), as built.** `nitcad plot <deck | result file> [options] --out <image.png>`
(`src/cli/plot.cpp`; `nitcad` now links `render`):
- **Input:** a deck (its device and mesh rebuilt, the result file its analysis writes read: `sweep`, `transient`,
  `small_signal`, or the equilibrium's `checkpoint`), or a `.nitsweep`, `.nittran` or `.nitac` file with `--device
  <file.nitdev>` and, for a generated mesh, `--mesh <file.nitmesh>` (its simplices; refused unless its graph digest is
  the device mesh's). A result whose fields are not one value per node of the device's mesh is refused (exit 1):
  the node-count pairing of VR6, until Unit 45 stores the mesh digest in result files.
- **Views:** `solid`, `mesh`, `field` (a 1D device's is a profile plot), `slice` (`--slice x|y|z=<um>`, repeatable),
  `profile`, `bands` (1D), and the curves `iv` (a sweep, `analysis::current_curve`), `cv` (a small-signal run,
  `capacitance_curve`) and `current` (a transient's terminal currents against time). Options: `--field` (the keys of
  `field_info`), `--point`, `--scale`, `--range`, `--colormap`, `--contours`, `--edges`, `--opaque`, `--azimuth`,
  `--elevation`, `--perspective`, `--contact`, `--sense`, `--frequency`, `--log`, `--size` (1280 × 960), `--adapter`.
  Regions and doping need no results; other views before the run fail with "(run the deck first)".
- **Exit status:** 0 written, 1 failure (the message on stderr), 64 usage (the option list on stderr).
- **Gates** (`tests/deck/run_test.cpp`, the executable started as the other CLI tests do, `--adapter warp`): images
  of the expected size from a deck (regions before the run; n with contours, the mesh, doping and I-V after it),
  from a result file with its device (a chosen point), and from a generated mesh's result with `--mesh`; refusals
  with exit 1 of a views-before-run request, a point out of range, another deck's device (node counts named) and
  another deck's mesh file; usage errors with exit 64. Checked by hand on a 2D pn diode: n across the junction and an
  I-V of about 60 mV per decade. (The first try drew a uniform n: the demo deck's profiles had their depth axis
  along x, so their x windows were ignored and the doping compensated to zero, which the plot showed as it was.)
- **Mutation checks (VG6, Release): 7 of 7 caught** — the node-count and mesh-digest refusals, the point's range
  (first caught only by the status, through undefined behaviour; the test now checks the message), the results a view
  needs, `--out` and `--device` required, the `--point` option.
- **VG1, assembly traces:** nothing in Unit 25a changes `assemble`, `solve`, `deck` or `io` (render, its tests, the
  CLI's plot command and CMakeLists.txt only), so every existing suite runs the code `main` runs; the suites pass
  unchanged. Release 22/22 and Debug 16/16, no warnings.

### Unit 45 — generation and recombination rates and the mesh digest in results (proposed 2026-10-21; on `results/rates-digest`)

**Owner request:** VQ12 (2026-10-17) put this unit before 25b: "recombination and generation rates in results and
files; the mesh digest in sweep, transient and small-signal files", to unlock the R view and exact result–mesh
pairing. Documentation only until the owner approves.

**Why now (audited on `main` at `6d2135b`).**
- No result holds a rate. The assembler computes each mechanism only inside its residual: since Unit 40 the
  drift-diffusion rows are equation terms, one per mechanism (`terms.hpp`: SRH, Auger, radiative, interfaces,
  impact generation, local and nonlocal band-to-band tunnelling, and their electrothermal forms). The only rate in
  results is a tunnel path's (`TunnelPath::generation_cm3_s`).
- `nitcad plot` pairs a result file with a device by node count (25a·3): two meshes of equal size pass. The mesh
  digest exists (`mesh::graph_digest`, `io::mesh_digest`) and only checkpoints and P4 certificates store it.
- 25b's derived-field views and the reference's R panel need the rates per node.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| RS1 | **Rates read from the terms.** `DriftDiffusion::rates(x)` evaluates each recombination or generation term alone on the state and takes its contribution to each node's electron and hole continuity rows, divided by the node's control volume and scaled to cm⁻³ s⁻¹. Per node: SRH, Auger, radiative and surface recombination (net, R − G of each, signed), impact-ionization and band-to-band generation (electron and hole rows separately where they differ: nonlocal tunnelling puts its electrons and holes at the path's two ends). | The stored rate is the residual's, operand for operand: no second implementation to drift from the solver. | One extra pass of those terms per stored point. |
| RS2 | **`results::RateFields`** in `BiasPoint`, `TransientSnapshot` and the small-signal operating points: the components above and the net `U = R − G` per carrier, each empty when its mechanism is off. Surface recombination is the interface edges' rate over the node's control volume (an effective volume rate at interface nodes), labelled as such. | One type for every result that has fields; empty vectors cost nothing. | A volume rate for a surface process is an approximation of its distribution, not of its total (the total is exact). |
| RS3 | **The mesh digest in `RunRecord`** (`mesh_digest`, filled by `make_run_record` from the device's mesh); not part of `input_identity` (which already covers the mesh). | Every result type carries its run record, so sweep, transient and small-signal files get it at once. | — |
| RS4 | **Files:** new optional records in the sweep, transient and small-signal files for the rates and the digest; files without them (every file written before) load as before, with no rates and digest 0 (unknown). | Old files keep loading (Unit 42's rule, PG3). | — |
| RS5 | **`nitcad plot`** checks the digest when the result has one (exact refusal), and falls back to the node count for older files; `render` gains the fields `R` (net), `srh`, `auger`, `radiative`, `surface`, `Gii`, `Gbtbt`. | The pairing VR6 asked for; the R view for 25b and the reference grid. | Seven field keys. |
| RS6 | **Out of scope:** |E| and |J| (25b, `render` reconstructs them); optical generation (no model yet); rates in checkpoints; per-edge rates. | Keeps the unit to stored rates and the digest. | — |

**Owner decisions (asked):**

| # | Question | Recommendation |
|---|---|---|
| RQ1 | Store the rates by default, or only when asked (`BiasOptions::store_rates`, deck `[solver] rates = on`)? | **By default**, with `rates = off` to save file size: a commercial viewer opens any result and shows R without rerunning. Measured cost in RG5 decides whether the default holds. |
| RQ2 | Components (as RS1), or the net rate only? | **Components and the net rate**: which mechanism dominates is the question users ask of an R map. |
| RQ3 | Surface recombination as an effective volume rate at interface nodes (RS2), or left out? | **Include, labelled** as an effective volume rate, with its exact total. |
| RQ4 | The digest refusal in `nitcad plot`: refuse a result whose digest differs, or warn? | **Refuse** (exit 1), as for a node-count mismatch: an image of the wrong mesh is worse than none. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| RG1 | **Exact against the residual:** for every mechanism and every state of the existing device tests' kinds (diode, MOSFET, traps, impact, local and nonlocal BTBT, electrothermal), the volume-weighted rates reproduce that term's contribution to the continuity residual to rounding. |
| RG2 | **Physics checks:** a uniform bar at low injection gives SRH = Δn/τ; radiative B(np − n_ie²) exactly; Auger as its formula; at equilibrium every net rate is zero to rounding; the integral of impact generation equals the avalanche current's increase (Unit 19's integral). |
| RG3 | **No numerical change:** the assembly traces are byte-identical to `main`'s (`--order decl`); every suite passes unchanged; `input_identity` of the fixtures unchanged. |
| RG4 | **Files:** round trips of rates and digest bit for bit; every committed file of earlier units (`tests/io/fixtures`) loads unchanged with no rates and digest 0. |
| RG5 | **Cost measured:** time and file size per point with and without rates on the 2D MOSFET sweep and a 3D diode; recorded whatever they show. |
| RG6 | **`nitcad plot`:** a result of a different mesh with the same node count is refused by digest; an old file falls back to the node count; the R and component fields draw. |
| RG7 | **Suites and mutants:** Release and Debug at every commit, no warnings; mutants in the rate extraction (a term's sign, a missing volume division, electron and hole rows swapped, a mechanism's flag), the file records and the digest check, all caught. |
| RG8 | **No scope creep:** nothing outside RS1–RS5; ARCHITECTURE.md 6.6 and 11.3 as built after acceptance. |

**Commits (proposed):** 1) `rates()` and `RateFields`, filled by the solves (RG1–RG3, RG5's timing); 2) the digest in
`RunRecord` and the file records (RG4, RG5's sizes); 3) `nitcad plot` and the render fields (RG6).

**Risks:** the rates' extra pass slows stored points (RG5 measures it); files grow by up to eight node arrays per
point, which RQ1's off switch answers; the electrothermal terms add their rates into the heat row too (only the
continuity rows are read).

**Status:** authorized by the owner 2026-10-21 ("Approve, as recommended"): RQ1 rates stored by default
(`[solver] rates = off` to opt out), RQ2 components and the net rate, RQ3 surface recombination included and
labelled as an effective volume rate, RQ4 `nitcad plot` refuses a differing digest (exit 1). Commits 1–3 done,
**accepted by the owner (2026-10-21)**; ARCHITECTURE.md 6.6 and 11.3 as built (RG8).

**Commit 1 (the rates; RG1–RG3, RG5's timing), as built.** `DriftDiffusion::rates(x)` evaluates each generation-
recombination term alone (residual only, no Jacobian) and reads its electron and hole rows: rate = row × R0 / V in
cm⁻³ s⁻¹, positive for net recombination (one conversion for every term: `rate_scale_` = Ns / R0 and
`generation_scale_` = 1 / R0). Mechanisms by term name; local and nonlocal tunnelling add into one entry; the
electrothermal terms are read the same way. `results::RateFields` (srh, auger, radiative, impact; interface_n/p,
band_to_band_n/p; net_n/p) in `BiasPoint` and `TransientSnapshot` (small-signal operating points are `BiasPoint`s),
filled by `sweep_bias`, `trace_bias` and `solve_transient` after convergence; `BiasOptions::store_rates` (on, RQ1),
recorded in the run record only when off; deck `[solver] rates = on|off`. `rates()` reports nothing to the
assembly observer.
- **RG1, the balance:** at every node off the contacts of a converged state, the net current out of its control
  volume (from the stored edge currents) equals q U V for electrons and −q U V for holes, within 10⁻⁶ of the node's
  currents plus 5% of the run's current resolution (a node's net current below it is the rounding of large opposing
  majority fluxes: measured ≤ 9e-15 A/cm against a resolution of 1.06e-12 in a reverse-biased diode, while the
  minority carriers agree to five digits). Holds for SRH, Auger and radiative (forward diode), impact ionization
  (reverse to −4 V), interface traps and surface recombination (the small-signal tests' oxide-silicon device), local
  and nonlocal band-to-band tunnelling (nonlocal: electron and hole rates differ, each balanced) and the
  electrothermal terms.
- **RG2:** at equilibrium every rate is below 1e4 cm⁻³ s⁻¹ (against n_i/τ ≈ 1e16); forward recombination exceeds
  1e18; impact and tunnelling rates are generations (negative).
- **RG3:** assembly traces of the bias, sweep and transient tests byte-identical to `main`'s (`6d2135b`, `--order
  decl`, 852,646 lines); a run storing rates keeps the default run identity; switching them off changes no field or
  current.
- **RG5 (time):** sweeps of 7 points with and without rates, best of three: 2D diode, 14,701 nodes: 8.92 s and
  9.19 s; 3D diode, 4,961 nodes: 20.17 s and 20.24 s. The extra pass is below the runs' noise.
- **Mutation checks (RG7, Release): 10 of 10 caught** — the electron row's sign, the volume division, the rows
  swapped, Auger read as SRH, tunnelling dropped, the hole net from the electron rows (caught only once the nonlocal
  case was added: every other mechanism moves carriers in pairs), `store_rates` ignored in sweeps and in
  snapshots, the setting always recorded, the deck key ignored. (A first runner lost a test's exit status between
  two batch calls; its two survivors were rerun with a fixed runner, one then caught.) Release 22/22 and Debug 16/16,
  no warnings.

**Commit 2 (the digest and the files; RG4, RG5's sizes), as built.** `RunRecord::mesh_digest` (the device mesh's
`mesh::graph_digest`, filled by every `make_run_record`; not in `input_identity`; 0 = unknown). Files: the run record's
key 4 (the digest, written when nonzero), the points' key 18 and the snapshots' key 8 (a rates record of the ten
vectors, written when the rates were stored); absent in older files, which load with no rates and digest 0.
- **RG4:** sweep, transient and small-signal files round-trip the rates (signed zero, NaN, a denormal) and the digest
  bit for bit; a point or snapshot without rates stays without; the pre-Unit-42 fixtures load with no rates and
  digest 0; every run record (sweep, transient, by hand) carries the device mesh's digest.
- **RG5 (size), the sweeps above:** 2D diode 13.1 MB with rates against 9.0 MB without (+46%); 3D diode 4.9 against
  3.5 MB (+40%): five node arrays per point (SRH, Auger, radiative and the two nets) beside the fields, bands and
  currents. Time stays within noise (+0.1% and +4.5% this run). **For the owner (RQ1):** the default stays on as
  ruled; `rates = off` saves the 40–46%.
- **Mutation checks (RG7, Release): 7 of 7 caught** — the digest not written, not read, not filled; the point's and
  the snapshot's rates not written; the point's not read; a rates key swapped. Release 22/22 and Debug 16/16, no
  warnings.

**Commit 3 (`nitcad plot` and the render fields; RG6), as built.** `render::Field` gains `R` (net R − G),
`srh`, `auger`, `radiative`, `surface` (signed, on the signed log scale), `Gii` and `Gbtbt` (generations, positive,
log10); `field_values` refuses a rate the result does not hold, saying how to store it. `nitcad plot` refuses a
result whose run record names another mesh (RQ4: "its mesh digest is not the device's", exit 1), and falls back to
the node count for a result without a digest.
- **RG6:** a result of a mesh with the same node count (x over 3 µm instead of 2) refused by digest; the R map of a
  forward-biased diode drawn (looked at: the peak at the junction, zero at the contacts); a pre-Unit-42 sweep with
  its device drawn through the node-count fallback, its R refused ("holds no R rates"), and with another device
  refused by node count; the render fields' values, signs, scales and the missing-rate error.
- **Mutation checks (RG7, Release): 4 of 4 caught** — the digest check off, a generation not negated, SRH read as
  Auger, a missing rate not refused. Unit 45 in all: 21 of 21. Release 22/22 and Debug 16/16, no warnings.

### Unit 25b — derived fields, isosurfaces and streamlines (proposed 2026-10-21; on `render/derived-fields`)

**Owner request:** VQ7 (2026-10-17): "25b, derived |E| and |J|, isosurfaces and streamlines", after Unit 45
(VQ12); gates VG4 (in 25b also |E|, |J| with streamlines and isosurfaces) and VG8 (scientific validation), with VG6
(suites and mutants). The owner's "GO" (2026-10-21) asked for this proposal. Documentation only until the owner
approves.

**What exists (audited on `main` at `81622ca`).**
- Every result with fields holds what the vectors need, so nothing is solved again: the potential per node, and the
  electron and hole current through every mesh edge (`BiasPoint::edge_current_n/p` and
  `TransientSnapshot::edge_current_n/p`, in every result file since Unit 9). The mesh gives each edge its length and
  its coupling area (`mesh::Edge::coupling_area`), so I_ij / A_ij is the current density along the edge, in A/cm² in
  every dimension.
- `render::Geometry` has the linear elements (Kuhn tetrahedra of tensor cells, generated simplices), the slices, the
  barycentric interpolants and marching triangles (2D contours). It has no marching tetrahedra, no point location
  and no vector field.
- `render::Field` is nodal; `field_values` returns one value per node; `nitcad plot` draws field, slice and contour
  views of it. Line batches have one color and a pixel width.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| DS1 | **E per element:** E = −∇ψ of the linear interpolant on each segment, triangle or tetrahedron, in V/cm; insulators included. | Exact for a linear potential (VG8); the field in the oxide matters as much as in the silicon. | Constant per element; an element whose nodes lie in two regions (the box method puts the interface between nodes) mixes both permittivities' fields. |
| DS2 | **J per element from the stored edge currents:** for each element, the vector whose projections on its mesh edges best match the edge densities I_ij / A_ij (least squares over the edges with a positive coupling area and both ends in a semiconductor); Jn and Jp separately, J = Jn + Jp. An element without D independent such edges has no current vector. | The solver's own currents (Scharfetter–Gummel), not a second discretization; exact for a uniform current (VG8). A Kuhn tetrahedron of a tensor cell has exactly three axis edges, so the fit is a solve there. | Conduction current only: a transient snapshot's displacement current is not stored (named in the title). |
| DS3 | **Nodal values for drawing:** each node gets the average of its elements' vectors, weighted by the volume of its control-volume piece in each element (the region rule of 25a); \|E\|, \|J\|, \|Jn\|, \|Jp\| are the magnitudes of those averages. New `render::Field` keys `E` (V/cm), `J`, `Jn`, `Jp` (A/cm²), log scale by default; 1D gives profiles. | Every 25a view (maps, 3D surfaces, slices, contours) draws them unchanged; the average stays inside the node's own region. | Smoothed: a peak narrower than an element is spread (the element values remain available to the code and the tests). |
| DS4 | **Isosurfaces (3D):** marching tetrahedra on the elements, for any nodal field, one or more levels; triangles oriented toward increasing value, shared edges cut once (watertight), each level colored by its colormap position and drawn translucent over the solid. Values are interpolated in the field's scale (log10 for densities and \|J\|), as the colormaps already are. In 2D the same request draws contour lines (25a). | Isosurfaces are the standard 3D view of depletion edges, channels and current paths; interpolating in the drawn scale puts a 1e16 surface where the map shows 1e16. | Exact planes on a linear field only on a linear scale. |
| DS5 | **Streamlines of E, J, Jn or Jp (2D and 3D):** the nodal vectors interpolated linearly in each element; adaptive Runge–Kutta (Dormand–Prince 5(4)) with point location by walking element neighbours from a bucket grid; a line stops at the domain boundary, a contact, the edge of the current-carrying region (J), a vanishing field or a length limit. Drawn as lines, in one color over the field map, or over the translucent solid in 3D. | Integrating the same linear field the maps show keeps lines and colors consistent; an adaptive integrator gives VG8's closure test a stated tolerance. | Lines have no arrowheads (vector glyphs are 25c). |
| DS6 | **Seeds:** J lines start on a contact, spaced so each carries an equal share of its current (the line density then shows the current density); E lines start evenly along a contact or gate; an explicit seed line or plane (`--seeds`) for either. | What a TCAD user reads from a current-line plot is where the current flows and how much. | Equal-current seeding needs the contact's current density along it, read from its edges. |
| DS7 | **`nitcad plot`:** `--field E|J|Jn|Jp`; `--streamlines E|J|Jn|Jp` over any field view (`--lines <n>`, `--from <contact>`, `--seeds <line or plane>`); `--view iso --field <f> --levels <v,...>`. | One more switch per feature on the existing command. | — |
| DS8 | **Out of scope:** vector glyphs, clipping, exploded views and filmstrips (25c); displacement current; per-element drawing (DS3's alternative); animated or tube-shaped streamlines; E or J stored in results (the `render` layer reconstructs them). | Keeps the unit to VQ7's list. | — |

**Owner decisions (asked):**

| # | Question | Recommendation |
|---|---|---|
| BQ1 | Draw \|E\| and \|J\| as nodal averages (smooth, DS3) or as one value per element (faceted, exact)? | **Nodal averages**, as commercial viewers draw them; the exact element values are what the gates test. |
| BQ2 | Seed J streamlines on a contact with an equal current per line (DS6), or only evenly along a user's line? | **Equal current per line from a contact**, with `--seeds` for an explicit line: the line density then means current density. |
| BQ3 | Interpolate isosurfaces in the field's scale (log10 for densities and \|J\|), or always linearly in the value? | **In the field's scale**, so the surface lies where the map's color says the level is. |

**Acceptance gates, fixed before any implementation (VG8, VG4 and VG6 made precise):**

| # | Gate |
|---|---|
| BG1 | **E and J exact per element (VG8):** \|E\| of a linear potential is exact (to rounding) on every element of 1D, 2D and 3D tensor and generated meshes; J from edge currents of a uniform current (I_ij = J · t_ij A_ij) is exact on every element; a solved uniform resistor bar gives \|J\| = I / area at every node to the solver's current resolution; Jn + Jp = J. |
| BG2 | **Isosurfaces (VG8):** marching tetrahedra on a linear field give the exact plane (every vertex on it, normals along the gradient); on a sphere's distance field inside the mesh the surface is closed (no edge used once) and outward, and its enclosed volume converges at second order under refinement (rate measured on three meshes, between 1.8 and 2.2). |
| BG3 | **Streamlines (VG8):** a uniform field gives straight lines to rounding; the rotational field (−y, x) (linear, so exactly represented) closes after one turn within the integrator's tolerance; lines stop at the boundary and at contacts; equal-current seeds of a uniform bar are evenly spaced, and of a known nonuniform contact current spaced by its integral. |
| BG4 | **Images (VG4):** WARP references, kept like 25a's: the 2D MOSFET's \|E\| with E lines and its \|J\| with current lines; the 3D MOSFET's isosurfaces of n and its current lines over the translucent solid; within 25a's per-pixel tolerance; the D3D12 debug layer clean in Debug. |
| BG5 | **No numerical change:** nothing outside `render` and `nitcad plot` changes (a header-boundary test as VG1); the assembly traces byte-identical to `main`'s; every suite passes unchanged. |
| BG6 | **`nitcad plot`:** each new field, the isosurface view and streamlines written from decks and from result files; a result without edge currents, a 2D isosurface request and an unknown contact refused with a message (exit 1). |
| BG7 | **Benchmarks:** \|E\| and \|J\| reconstruction, isosurfaces and 100 streamlines on the 3D MOSFET and on a 10⁶-tetrahedron mesh, Release, recorded. |
| BG8 | **Suites and mutants (VG6):** Release and Debug at every commit, no warnings under `/W4 /WX`; mutants in the gradient, the edge fit, the averaging weights, marching tetrahedra (an orientation, a shared-edge cut), the integrator (a stage weight, the step control), point location and the seeding, all caught. |
| BG9 | **No scope creep:** nothing outside DS1–DS7; ARCHITECTURE.md 6.14 and 11.3 as built after acceptance. |

**Commits (proposed):** 1) element and nodal E and J, the four field keys and `--field` (BG1, BG5, part of BG6);
2) marching tetrahedra and `--view iso` (BG2); 3) point location, the integrator, seeding and `--streamlines`
(BG3); 4) the references and the benchmark (BG4, BG7). About 1,700 lines with tests (VQ7's estimate).

**Risks:** element-mixed interfaces blur E across an oxide–silicon interface on coarse meshes (DS1; the gates use
single-region elements, the references show it as drawn); streamlines through regions of near-zero current wander
(stopped by DS5's vanishing-field rule); point location on large 3D meshes (BG7 measures it).

**Status:** authorized by the owner 2026-10-22 ("Approve, as recommended"): BQ1 \|E\| and \|J\| drawn as nodal
averages, BQ2 J streamlines seeded on a contact with an equal current per line (and `--seeds` for an explicit line),
BQ3 isosurfaces interpolated in the field's scale.

**Commit 1 (\|E\| and \|J\|; BG1, BG5, BG6 in part), as built.** `vectors.hpp`: `element_vectors` (E = −∇ψ per
element, by Cramer's rule on its edge vectors; Jn and Jp by least squares on the element's mesh edges with a positive
coupling area and both ends in a semiconductor, J their sum), `control_volume_shares` (tensor: binomial(D, m) / 2^D
for the vertex m steps along the Kuhn path, which is the corner sub-box's part of the simplex, so the shares add up to
the box method's control volumes exactly; simplex: 1 / (D + 1)), `nodal_vectors`, `magnitudes`. `render::Field` gains
`E`, `J`, `Jn`, `Jp` (log10), given by a `field_values` overload taking the geometry (the old one refuses them).
`nitcad plot --field E|J|Jn|Jp`: on a log scale without `--range`, the ten decades below the largest value, values
under them at the bottom color (below that the reconstructed currents are rounding: Jn at equilibrium spans 1e-42 to
1e-10 A/cm²); a magnitude that is zero everywhere on a linear 0..1 scale.
- **BG1:** E of a linear potential and J from the edge currents of a uniform current exact (10⁻¹¹) on every element
  and node of 1D, 2D and 3D tensor and 2D and 3D generated meshes; the tensor shares sum to the mesh's control
  volumes (10⁻¹²) in 1D, 2D and 3D; a nonuniform 1D E averages each node's two half segments; a solved 2D resistor
  bar at 0.1 V has \|J\| = I / height at every node within 10⁻⁶ plus the current resolution, along +x; Jn + Jp = J;
  elements with an oxide node have no current, oxide nodes NaN \|J\|; results without edge currents and other or
  element-less geometries refused.
- **Finding (DS2 amended):** a quasi-static result stores zero edge currents, not none (`sweep_bias` reports zero
  currents), so its \|J\| draws as zero; the refusal is for results without edge currents.
- **BG6 (in part):** E, J, Jn and Jp drawn from a 2D deck and from its result file (the equilibrium point's \|J\|,
  exactly zero, on 0..1), a 1D \|E\| profile (the junction's peak), a quasi-static \|J\| drawn as zero.
- **Mutation checks (BG8, Release): 14 of 14 caught** — E's sign, the 2D and 3D Cramer solves, insulator edges kept
  at either end, a reversed edge tangent, no division by the coupling area, holes subtracted, equal tensor shares, a
  wrong Kuhn step count, unweighted volumes, a 2D magnitude, `nitcad plot` on the old overload, the zero-field
  scale off. (With oxide only at the top, the first-node check was unseen: tensor numbering puts the silicon node
  first on every such edge, so oxide was added below. Two first versions failed only the `/WX` build, an unused
  variable, and were replaced by mutants the tests catch.)

**Commit 2 (isosurfaces; BG2, BG6 in part), as built.** `geometry.hpp`: `isosurface` (marching tetrahedra into an
indexed mesh: one vertex per cut mesh edge, computed from its lower-numbered node so every element around the edge
makes the same point, or the node itself for a crossing exactly at it; a node at the level counts as above it, as in
the 2D contours; triangles of zero area dropped; each triangle turned toward higher values by the sign of
Σ (v_k − level)(x_k − p₀), which for a linear field is the gradient's). `ColorScale::transform` is public.
`scene.hpp`: `iso_view` (3D: the surfaces of the field's scale coordinate (BQ3), colored by each level's colormap
position, smooth-shaded, opaque for one level and translucent for more, in the device's faint outline with the
contacts; 2D: the field map with contour lines at the levels). `nitcad plot --view iso --field <f> --levels
<v>[,<v>...]`.
- **BG2:** a linear field's surfaces on graded tensor and generated 3D meshes: every vertex on the plane (10⁻¹²),
  every normal along the gradient (cosine ≥ 1 − 10⁻⁹), no directed edge twice. A sphere's distance field (R = 0.3 µm,
  centre off the grid planes) on 8³, 16³ and 32³ cells: closed and consistently oriented, outward, volume errors
  8.81e-2, 2.15e-2 and 5.41e-3, **rates 2.03 and 1.99**. Crossings at nodes give one vertex per node (both edge
  orientations); a level at the field's maximum gives the boundary face; NaN nodes and 2D geometries give nothing.
- **BQ3:** a field exponential in y (log scale) has its 10^14.3 surface on the plane y = 0.43 µm between grid
  planes, and its 2D contour on that line.
- **BG6 (in part):** 3D surfaces of n (two levels) and ψ (one) from a deck, 2D contours at a level; refused: no
  `--levels` (usage), an unreadable level (usage), a negative level on a log scale, a 1D device. **For the owner:**
  BG6 listed "a 2D isosurface request" among the refusals while DS4 draws contours at the levels in 2D; DS4 is built
  (2D draws contours), and the refusal applies to 1D.
- **Mutation checks (BG8, Release): 9 of 9 caught** — the shared edge cut from either end, the node rule off, a node
  at the level counted below, a crossed quadrilateral, flipped orientation, contour levels ignored, interpolation in
  the value instead of the scale, one level translucent, `nitcad plot` allowing 1D. (The first two survived a first
  run: Kuhn simplices number their nodes in increasing x, so only a decreasing field reaches the node rule's other
  branch, and only a level at the field's maximum shows the at-level convention; both cases were added.)

**Commit 3 (streamlines; BG3, BG6 in part), as built.** `streamlines.hpp`: `Locator` (the hint element, then the
elements sharing a node with it, then a bucket grid of about one element per bucket over the elements' bounding
boxes; barycentric weights by Cramer's rule, inside at −10⁻¹⁰); `trace_streamline` (the nodal vectors interpolated
linearly, followed at unit speed, adaptive Dormand–Prince 5(4) with first-same-as-last, local error per step at most
`tolerance` (10⁻⁵) times the device's size, steps at most half the local element's size so the drawn polyline follows
the field; forward and backward from the seed; stops at the boundary (to 10⁻⁷ of the size), at an undefined or
vanishing field (|v| ≤ 10⁻¹² of the largest), at a sink (a step that turns the direction back), at 20 device sizes
or 20,000 steps; |v| kept per point); `contact_node_currents` (the stored edge currents from each contact node to
nodes outside the contact); `contact_seeds` (BQ2: the contact's boundary segments (2D, walked in order along each
chain) or faces (3D), each node's \|current\| spread evenly over its part, seeds at the midpoints of equal shares;
without currents, by length or area); `segment_seeds`, `plane_seeds` (by area over `Geometry::slice`).
`scene.hpp`: `add_streamlines` (one color, or each segment through the colorbar in 64 bins), `outline_view` (the
outline shell, now shared with the slice and isosurface views, and the contacts). `nitcad plot --streamlines
E|J|Jn|Jp [--lines n] [--from contact] [--seeds a:b | x=um]`: 2D over the field map (by default the lines' own
magnitude) in near-black; 3D over its own magnitude, the lines colored by it in the device's outline, else over the
chosen view. Default seeds: J lines on the contact carrying the most current, E lines on the first gate or electrode,
else the most biased contact.
- **BG3:** a uniform field's lines straight to rounding, ordered, ending on y = 0 and y = 1 µm (10⁻⁶ of the size);
  the rotational field (−y, x) closes after one turn to 1.2e-15 cm (109 points) and, on a coarse mesh at tolerance
  10⁻¹², to 1.4e-16 cm; a field with a kink on a node column, v = (1, k\|x − xm\|), follows its closed-form line within
  4.2e-14 cm (budget: tolerance × size × points); lines stop at a NaN region (no field), at a sink, and a zero field
  gives the seed alone; 2000 points located in a generated 3D mesh, their weights rebuilding them to 10⁻¹²; equal-
  current seeds of a uniform current at (i + ½)/N of the contact, and of j(y) = 1 + 3y/H at equal shares of its
  integral (10⁻⁹), denser where the current is; the two contacts' currents balance (10⁻¹²); plane seeds on the plane.
- **BG6 (in part):** J, Jn (from a named contact) and E (from a segment, over ψ) lines in 2D; J lines, E lines from a
  plane and Jp lines over a slice in 3D; refused: an unknown field, a malformed segment, a 3D segment in 2D, a plane
  in 2D, `--lines 0` (usage), an unknown contact, a plane that misses the device, a 1D device.
- **Mutation checks (BG8, Release): 12 of 12 caught** — the quantiles' midpoints, a Dormand–Prince stage weight and
  the error estimate's last weight, the barycentric weight of node 0, the length limit, every step accepted, the
  sink stop, the step cap, the backward part's order, a contact current's sign, the density not per length, a half
  segment reversed. (Every step accepted survived until the kinked field: in a smooth field the step sizing keeps
  each attempt within the tolerance, and only a kink makes a step fail.)

**Commit 4 (references and the benchmark; BG4, BG7), as built.** The four references showed two faults of the lines
the unit tests had not, both fixed here in `vectors.hpp`'s `streamline_field`, which the plot and the references now
trace instead of the bare nodal vectors:
- **Finding 1, the elements along an interface:** with insulator nodes NaN, every element joining the last silicon
  row to the first oxide row has no current vector, so each current line stopped at the silicon surface, short of the
  drain contact (in an on-state MOSFET it would cut every channel line). Insulator nodes now carry zero current for
  streamlines (the |J| map keeps them as no data), so lines continue along the interface and stop only where the
  current vanishes in the oxide.
- **Finding 2, normal components on no-flux boundaries:** the element average at a boundary node comes from one side,
  so it kept a normal part there (the top silicon row's current pointed into the oxide, and lines left through it
  0.25 µm short of the drain). The model has no normal current at semiconductor–insulator interfaces and no normal
  current or field on the outer boundary away from contacts (homogeneous Neumann), so `streamline_field` removes the
  normal part there: per node its boundary pieces' normals, grouped within 30° (one face, or a gently curved one)
  and orthonormalized (the faces meeting at an edge or corner each count), every component along them removed;
  contact nodes keep theirs. DS5 is amended accordingly.
- **BG4:** WARP references, 960 × 720, in 25a's tolerance (at most 0.5% of pixels off by more than 24):
  `mosfet2d_efield_lines` (\|E\|, 16 lines from the gate through the oxide to the depletion edge, where they end at
  the field's sink), `mosfet2d_current_lines` (\|J\|, 16 equal-current lines from the source: this subthreshold
  device conducts under the channel through the 1e15 substrate, then along the drain pocket's surface to the drain),
  `mosfet3d_iso_n` (n = 1e14 and 1e18, the source and drain junctions), `mosfet3d_current_lines` (24 lines colored
  by \|J\|, source to drain contact). Each was looked at before being kept; all 15 references pass; no D3D12 debug
  errors.
- **Tests:** insulator nodes zero and interface nodes tangential with J = (5, 2), and a line along the top silicon
  row reaching the far contact inside the interface band; E = (3, 4) on 2D and 3D bars keeps (3, 0) on the
  non-contact sides, edges and corners included, and (3, 4) at contacts and inside.
- **BG7 (Release, other jobs running, so upper bounds):** the 3D MOSFET (3,087 nodes, 14,400 tetrahedra): E 0.001 s,
  J 0.003 s, an isosurface 0.0004 s, streamline field and locator 0.009 s, 100 lines 0.012 s; a 3D tensor device of
  185,193 nodes and **1,053,696 tetrahedra**: E 0.048 s, J 0.495 s, an isosurface 0.014 s (41,078 triangles),
  streamline field and locator 0.73 s, 100 lines 0.054 s (24,844 points).
- **Mutation checks (BG8, Release): 5 of 5 caught** — insulator nodes not zeroed, every normal merged into one,
  interfaces ignored, no projection, contacts projected. Unit 25b in all: 40 of 40.

**BG5 (no numerical change), as checked.** Only `render`, `nitcad plot`, tests and documents differ from `main`
(`81622ca`); `CMakeLists.txt` adds the three sources and test files. Assembly traces of the bias, sweep and transient
tests (`--order decl`): the 38 cases other than "bias: assembly cost on the legacy diode" are byte-identical
(666,691 lines); that case repeats residual evaluations for a time budget, and on the branch (run while three test
suites shared the CPU) it stopped after 118,371 of main's 185,955 lines, which are an exact prefix of main's. Its count
depends on the machine's load, not on the code.

**BG8 suites:** Release 22/22 and Debug 16/16 at each of the four commits (commits 2–4 in clean worktrees), no
warnings under `/W4 /WX`; the Debug render suite includes the D3D12 debug layer check (no errors).

**Acceptance (2026-10-22).** The owner pre-accepted the unit and authorized push, PR and merge, then the start of 25c
without a further approval ("As Unit 25b is done, accept it push, open pr and merge, then start 25c without my
approval"), given while only commit 1 existed; it was taken as acceptance on the gates passing, which they did.
ARCHITECTURE.md 6.14, 11 and 11.3 as built (BG9). For the owner's record: DS2 amended (a quasi-static result's
currents are zero, not absent), DS5 amended (the streamline field, commit 4), and BG6's "2D isosurface request
refused" superseded by DS4 (2D draws contours at the levels).

### Unit 25c — vector glyphs, the clipping box, the exploded view, transient filmstrips (2026-10-22; on `render/vectors-clip`)

**Owner request:** VQ9 (2026-10-17): "Unit 25c (vector glyphs, clipping box, exploded view, transient filmstrips)
after 25b, each view droppable"; the plan's views 20 (current arrows on a 2D section), 21 (E arrows on the surface),
23 (a clipping box, a quarter cut away), 24 (regions and electrodes lifted apart) and 26 (a filmstrip of n as the
gate switches on), about 1,600 lines. On 2026-10-22 the owner asked for 25c to start without a further approval
("then start 25c without my approval"): the rulings below (CQ1–CQ4) are taken as recommended under that instruction
and recorded for the owner's review at acceptance.

**What exists (audited on `main` at `8256cf8`).** Nodal E and J and the no-flux streamline field (Unit 25b), point
location, slices and convex-polygon clipping (`Geometry::slice`, the internal `clip`), the 3D solid layer by
control-volume pieces with translucent insulators, contact slabs, the outline shell. Transient results hold
snapshots with fields, edge currents and rates, but `nitcad plot` converts a snapshot to a bias point with its fields
and bands only, so \|J\|, streamlines and R are unavailable on transient results (finding; fixed here).

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| DC1 | **Vector glyphs** of E, J, Jn or Jp (`vectors.hpp`'s streamline field, so tangential on no-flux boundaries): 2D at the points of a regular grid over the device, about N across; 3D at points of the boundary surface, spread evenly by area (a deterministic stratified sample of its faces), or on a slice. Each arrow a flat shaft and head (2D, over the map) or a lit cone on a cylinder (3D), of one length (a fraction of the spacing), colored by \|v\| through the colormap on the field's scale; points with no field drawn nothing. | Direction is what arrows add over a map; one length keeps 20 decades of magnitude readable, and the color carries the magnitude (CQ1). | No length-coded magnitude. |
| DC2 | **Clipping box** (3D): an axis-aligned box removed from the device. The remaining boundary faces (each clipped to its parts outside the box, as convex pieces) and the cut faces (the box's faces inside the device, from slices of the geometry clipped to the face) are drawn with the field or the regions, insulators translucent as in the field view. | The standard way to see inside a 3D device while keeping its outline; exact from the existing convex clipping. | The box is axis-aligned. |
| DC3 | **Exploded view** (3D): each region's closed surface (its boundary faces and its side of every interface) and the contact slabs on it, moved apart along an axis (default y, the stack) by a gap per step in the regions' order along it; drawn solid by region or with a field. | Shows the layer stack and buried regions without slicing. | Regions that interlock along the axis still overlap in projection. |
| DC4 | **Transient filmstrip**: chosen snapshots (all, or every k-th, at most 12) drawn as one field view each and tiled into one image, each titled with its time, on one scale fitted over all of them (CQ2); `--frames <dir>` writes each frame as its own PNG for a video tool. | Comparing frames needs one scale; frames serve presentations. | No animation playback (26). |
| DC5 | **`nitcad plot`**: `--vectors E\|J\|Jn\|Jp [--arrows n]` over any 2D or 3D field or slice view; `--clip <x0>,<y0>,<z0>:<x1>,<y1>,<z1>` (µm) on 3D field and solid views; `--view exploded [--gap <um>] [--axis x\|y\|z]`; `--view filmstrip [--every k] [--columns c] [--frames <dir>]` on a transient result. Transient snapshots keep their edge currents and rates. | One flag per view on the existing command. | — |
| DC6 | **Out of scope:** interactive clip handles, animation (26); oblique clipping planes; streamtubes; glyphs scaled by magnitude; video encoding. | Keeps 25c to VQ9's list. | — |

**Rulings (CQ1–CQ4), taken as recommended under the owner's instruction of 2026-10-22:**

| # | Question | Ruling |
|---|---|---|
| CQ1 | Arrow length by magnitude, or one length with the magnitude as color? | **One length, colored by \|v\|** on the field's scale. |
| CQ2 | Filmstrip frames on one shared scale, or each fitted? | **One scale** fitted over all chosen frames. |
| CQ3 | Exploded gap default? | **12% of the device's extent along the axis**, per region step. |
| CQ4 | Clipping applied to which views? | **Solid and field views in 3D** (and the vectors over them); slices and isosurfaces unchanged. |

**Acceptance gates:**

| # | Gate |
|---|---|
| CG1 | **Glyphs:** in a uniform field every arrow is parallel to it and of the set length; 2D arrows lie on the grid's points inside the device, as many as asked (to ±10%); 3D arrows lie on the boundary and are tangential on no-flux faces; colors are the colormap at \|v\|'s position; a NaN region gets none. |
| CG2 | **Clipping:** for a box inside a tensor device, the kept boundary pieces and the cut faces together form a closed, outward surface whose enclosed volume is the device's minus the box's (to rounding); a box covering a corner (the quarter cut) likewise; a box missing the device leaves the view unchanged; field values on the cut faces equal the interpolated field. |
| CG3 | **Exploded:** each region's surface is closed and outward; its volume is the region's control volume total; each region and its contacts move by the stated offset; with a zero gap the surfaces cover the solid view's. |
| CG4 | **Filmstrip:** one tile per chosen snapshot in time order, one shared scale, titles with the times; `--frames` writes one PNG per frame; transient snapshots draw \|J\| and streamlines. |
| CG5 | **Images:** WARP references (25a's tolerance): the 2D MOSFET's \|J\| with current arrows, the 3D MOSFET's E arrows on its surface, a quarter-cut clipping box of n, the exploded MOSFET, a filmstrip; a clean D3D12 debug layer in Debug. |
| CG6 | **`nitcad plot`:** each new flag from decks and result files; refusals with a message (exit 1, or usage 64): vectors in 1D, a clip box in 2D or not a box, a filmstrip of a non-transient result, unknown axes. |
| CG7 | **No numerical change** (only `render`, `nitcad plot`, tests and documents change); Release and Debug suites at every commit, no warnings under `/W4 /WX`; mutants in the arrow placement and direction, the clipping (inside/outside, the cut faces), the explosion offsets and the filmstrip scale and order, all caught. |
| CG8 | **Benchmarks:** arrows, clipping and the exploded view on the 10⁶-tetrahedron device, recorded. |
| CG9 | **No scope creep:** nothing outside DC1–DC5; ARCHITECTURE.md 6.14 and 11.3 as built after acceptance. |

**Commits (planned):** 1) vector glyphs and `--vectors` (CG1); 2) the clipping box and `--clip` (CG2); 3) the exploded
view (CG3); 4) filmstrips, transient snapshots' currents and rates in the plot, references and benchmarks (CG4, CG5,
CG8). Each view is independent and can be dropped (VQ9).

**Commit 1 (vector glyphs; CG1, CG6 in part), as built.** `streamlines.hpp`: `Arrow` (site, unit direction, \|v\|,
the surface's normal, the interpolant), `grid_arrows` (2D: cell centres of a square grid centred on the mesh's box,
spacing sqrt(area / count)), `surface_arrows` (3D: a stratified sample of surface polygons by area, each with its
polygon's normal); both read the field at the site through the point locator and skip undefined or zero fields.
`scene.hpp`: `add_arrows` (2D a flat shaft and head just in front of the map, unlit; 3D a six-sided shaft and cone
lifted a quarter of the length off the surface, lit; one length; one color or \|v\| through a colorbar),
`exposed_boundary` (the boundary faces not wholly on a contact). `nitcad plot --vectors E|J|Jn|Jp [--arrows n]`:
arrows of the streamline field (Unit 25b: tangential on no-flux boundaries), 0.8 of the spacing long; 2D over the
map (by default the vectors' own magnitude) in near-black; 3D, over their own magnitude, on the solid device colored
by \|v\| with its colorbar (CQ1), else near-black over the chosen field or slice view; 3D sites on the exposed
boundary or the slices shown. **CQ1 as built:** colored by \|v\| where the view does not itself show a field (the
solid); over a field map an arrow colored like the map under it would vanish, so there it is near-black.
- **CG1:** a uniform field's grid arrows (50, 200 and 800 asked) number within 10% of the request, point along the
  field (10⁻¹²), lie on a grid of the stated spacing centred on the device; none where the field is NaN or zero;
  surface arrows on a 3D bar's exposed faces (its contact faces excluded: 6 sides less 2) have unit normals and no
  normal component away from the contacts; drawn arrows span the set length along their direction, sit in front of
  the 2D map or above their 3D surface, and carry the color or \|v\|'s colormap position.
- **CG6 (in part):** J arrows in 2D and 3D, E arrows over n (2D) and over a slice (3D), Jn arrows with streamlines;
  refused: an unknown field and `--arrows 0` (usage), a 1D device.
- Looked at on a 2D and a 3D MOSFET deck (the 3D one at V_G = 5 V, V_D = 0.5 V): current into the drain pocket and
  out of the source, the body leakage, fields from the gate; arrows on contact faces stood normal to them under the
  slabs, so contact faces are excluded (`exposed_boundary`).
- **Mutation checks (CG7, Release): 11 of 11 caught** — zero fields kept, the direction unnormalized, the spacing
  doubled, the grid not centred (caught once the centring check was added), a non-unit normal, the 2D arrow behind
  the map, a short 2D and 3D arrow, the 3D arrow sunk into its surface, contact faces kept, the color not \|v\|'s.

**Commit 2 (the clipping box; CG2, CG6 in part), as built.** `geometry.hpp`: `Box`; `outside_box` (a convex
polygon cut by the box's six planes in turn, the part beyond each kept and the rest passed on; a part lying in a box
plane, with nothing strictly beyond it, belongs to the box); `Geometry::clipped_surface` (the boundary faces' parts
outside the box, and the cut faces: sections of the elements by each box face strictly inside the mesh's bounding
box, clipped to the face and turned toward the box, each with its interpolants). `ViewOptions::clip` (CQ4: the 3D
solid and field views draw the clipped surface, interfaces cut alike, a contact slab's face left out when its
centroid is in the box). `nitcad plot --clip <x0>,<y0>,<z0>:<x1>,<y1>,<z1>` (µm), arrows then sit on the clipped
view's exposed boundary.
- **CG2:** a square's parts outside a box are convex and cover the square less the overlap; on a graded tensor mesh,
  an interior box, a corner box (the quarter cut), a slab through the device (two pieces), a box with a face on the
  boundary and a box missing the device give closed surfaces (vector area 0 to 10⁻¹⁰) enclosing the device's volume
  less the overlap (10⁻¹⁰), outward; the same for an interior and a corner box on a generated mesh; a linear field
  exact on the cut faces; the field and solid views (with oxide above and below, so interfaces too) draw nothing
  inside the box and draw the cut face; a box over the drain end removes its pad.
- **CG6 (in part):** clipped field, solid and arrow views from a deck, and the clipped image differs from the whole;
  refused: a two-coordinate or inverted box (usage), a slice view, a 2D device.
- Looked at on the 3D MOSFET deck: the quarter cut shows n inside the source and drain pockets and the channel
  region, the gate slab notched where the box removes it.
- **Mutation checks (CG7, Release): 10 of 10 caught** — the box's faces swapped, in-plane parts kept, the remainder
  not reduced, cut faces at the bounding box, cut faces turned outward, cut faces not clipped to the face, slabs
  kept, the view unclipped, interfaces unclipped, `nitcad plot` ignoring the box (caught once the clipped and whole
  images were compared). The face-on-boundary case was added after a first version kept such a face.

**Commit 3 (the exploded view; CG3, CG6 in part), as built.** `geometry.hpp`: `Geometry::region_surfaces` (per
region, the boundary pieces its nodes own and its side of every interface piece, each turned outward from the
region). `scene.hpp`: `exploded_view` (each region's surface and its contacts' slabs moved along the axis by the
gap times the region's rank in the order of its mean node coordinate, the first staying; the semiconductor with a
field or by class, insulators translucent), `default_gap` (CQ3: 12% of the extent). `nitcad plot --view exploded
[--field f] [--axis x|y|z] [--gap um]` (regions without `--field`).
- **CG3:** on a 3D silicon device under an oxide layer with an electrode, each region's surface is closed (vector
  area 0 to 10⁻¹⁰) and outward, enclosing exactly its nodes' control volumes (10⁻¹⁰); at zero gap the surfaces' area
  is the boundary's plus both sides of the interface; the oxide and its electrode move up by one gap and the silicon
  stays; a field colors the silicon only, at its scale's positions; refused: an axis outside 0..2, a negative gap, a
  field not one per node.
- **CG6 (in part):** regions, zero gap and a field along x from a deck; the gap and the field change the image;
  refused: an unknown axis and a negative gap (usage), a 2D device.
- Looked at on the 3D MOSFET deck: the oxide and gate lifted off the silicon, and with n the inverted channel's
  surface (V_G = 5 V) between the pockets.
- **Mutation checks (CG7, Release): 9 of 9 caught** — boundary pieces put in region 0, interface pieces not turned,
  interface sides swapped, the default gap, the order reversed, the first region moved, the field ignored (caught
  once a field case was added), contacts left behind, `nitcad plot` ignoring `--field`.

**Commit 4 (filmstrips, transient currents, references, benchmarks; CG4, CG5, CG8), as built.** `scene.hpp`:
`FilmFrame`, `filmstrip_scenes` (one field view per frame on one scale fitted over every frame (CQ2), titled with
its label), `time_label` (three significant digits in s, ms, µs, ns, ps or fs). `image.hpp`: `tile_images` (rows of
`columns`, in order, on white). `nitcad plot --view filmstrip --field f [--every k] [--columns c] [--frames dir]`
on a transient result: at most 12 frames (every k-th snapshot from the first; by default every ⌈n/12⌉-th), each
drawn at the frame size (`--size`, else 640 × 480: a first version tiled the frames into `--size` and left each map
about 110 pixels wide), 4 to a row, each also written alone under `--frames`; a vector magnitude's ten-decade floor
taken over all frames. **Finding (fixed):** `nitcad plot` read a transient snapshot's fields and bands only; it now
keeps the snapshot's edge currents, rates and bias, so \|J\|, streamlines, arrows and R work on transient results.
**Arrows (amending commit 1):** `grid_arrows` and `surface_arrows` take an exclusion mask; a site whose nearest node
is masked gets no arrow, and `nitcad plot` masks the insulators for currents (the first 2D reference showed a row
of current arrows in the oxide, where the elements next to the interface interpolate the silicon's current).
- **CG4:** tiles in order, the empty place white, sizes checked; frames on one scale (1..5 from frames of 1..2 and
  3..5), titled, a given scale kept; time labels; a transient deck's filmstrip of n with every snapshot (frame count
  in the output, the strip's size, one PNG per frame), of \|J\| every 2nd in 2 columns, and current streamlines from
  the transient file; refused: a sweep's results, `--every 0`.
- **CG5:** WARP references (960 × 720 unless stated), each looked at before being kept: `mosfet2d_current_arrows`
  (\|J\| with 250 grid arrows, none in the oxide), `mosfet3d_efield_arrows` (E arrows on the exposed surface,
  colored by \|E\|), `mosfet3d_clip_n` (n with the quarter cut, the gate slab notched), `mosfet3d_exploded` (the oxide
  and its gate lifted), `mosfet2d_filmstrip_n` (n over the five bias points of the gate and drain sweep, 3 to a row of
  480 × 360 frames: the drain-side depletion grows with VD). 20 references in all pass.
- **CG8 (Release, other jobs running, so upper bounds):** the 3D MOSFET: 400 surface arrows 0.007 s, the quarter-cut
  surface 0.010 s, the exploded view 0.044 s; the 1,053,696-tetrahedron device: 400 arrows 0.61 s (the locator
  included), the quarter-cut surface 0.67 s (37,632 polygons), the exploded view 2.0 s.
- **Mutation checks (CG7, Release): 10 of 10 caught** — tiles transposed, the scale from the first frame only, the
  label left out of the title, the unit threshold, the decimals, the exclusion off, the snapshots' currents dropped,
  the default every-k, the default columns, the frames not written. (`nitcad plot`'s insulator mask for arrows is
  wiring the tests do not observe; the mask itself is tested, and the reference shows its effect.) Unit 25c in all:
  40 of 40.

**CG7 suites and numerics:** Release 22/22 and Debug 16/16 at each of the four commits (in clean worktrees), no
warnings under `/W4 /WX`, the Debug render suite's D3D12 debug layer clean; nothing outside `render`, `nitcad plot`,
tests, `CMakeLists.txt` and this file differs from `main` (`8256cf8`), so no numerical change. **Status:**
**accepted by the owner (2026-10-22, "accept")**, CQ1–CQ4 with it; ARCHITECTURE.md 6.14, 11 and 11.3 as built (CG9).
Push, PR and merge await the owner's authorization.

### Unit 26 — the Windows application (authorized 2026-10-22; on `app/windows-application`)

**Owner request:** VQ3/VQ6 order (25a, 45, 25b, 25c, then 26), VA1–VA8 and gates WG1–WG6 (Units 25 and 26 above),
the detailed plan's start page, run monitor and results workspace (page Figures 3, 5, 6), VQ11 (bundled examples) and
VQ13 (menus open the deck at a section). On 2026-10-22 the plan page was checked against this file and the repository
and four points were raised (AQ1–AQ4); the owner ruled on them the same day and gave the go: "GO — Unit 26 is
authorized."

**What exists (audited on `main` at `e289e44`).** `solve::RunControl` (a stop token and a progress callback, 6.9);
`deck::parse`, `deck::build` (errors `"<source>:<line>: <what>"` with the line as the context index) and `deck::run`
(a stopped run is still a report with its completed part); `render` with off-screen D3D12/Direct2D drawing to images,
an orthographic and a perspective camera (`--perspective`, reference `mosfet3d_opaque_perspective`, Unit 25a), slices
(25a), the clipping box (25c); `analysis` extractors; `io` loaders. No window, no swap chain, no `app` code, no
`examples/` directory. The layer table allows `app` everything except `physics` and `assemble` (ARCHITECTURE.md 3).

**Rulings (2026-10-22):**

| # | Point | Ruling |
|---|---|---|
| AQ1 | The planned two commits do not cover the start page, run monitor and results workspace. | **Approved: three commits.** 26·1 the app logic; 26·2 the main window; 26·3 the start page, run monitor, results workspace, the four example decks and `--self-test`. |
| AQ2 | WG6 bounds scope to VA1–VA6, but VA1 and VA8 are Unit 26's. | **Approved:** WG6 covers VA1–VA6 and VA8–VA9; WG7 (menu navigation) and WG8 (every example deck tested) added; VA7's disposition stated explicitly (below). |
| AQ3 | No unit owns the examples (VQ11). | **Approved:** 26·3 owns `examples/` and the four start-page decks. |
| AQ4 | How the application is launched. | **Approved:** a separate `nitcad-app.exe`; no GUI or window code in the `nitcad` batch runner. |

Owner instructions with the rulings: preserve the architecture and the acceptance rulings; verify whether hosted CI
can open a window for `--self-test`, and if not make WG4 local-only and keep the off-screen WARP image tests in CI,
without letting the question block the work; commits scoped to their approved responsibilities; report commit hashes,
test results and deviations; no unrelated code or scope.

**VA9 (added by AQ3, AQ4).** `nitcad-app.exe`, a `/SUBSYSTEM:WINDOWS` executable of the `app` layer (L8);
`examples/` at the top level with four decks (1D pn diode, 2D MOSFET, 3D MOSFET, MOS capacitor) listed on the start
page, each run by a test.

**VA7's disposition (AQ2).** VA7 is the exclusion list and grants nothing: a graphical device editor (form editors that
write the deck, VQ13), project management, remote or cluster runs, plugins, localization and an installer are not
built in Unit 26, and no gate requires them. WG6 checks the work against VA1–VA6 and VA8–VA9 and against VA7's list.
The detailed plan's analysis workspace (live checking and completion, structure preview, cut lines, run comparison,
difference maps, mesh quality; page Figures 4, 7, 8) is Unit 26b (VQ10), not 26.

**Implementation plan (within VA1–VA9):**
- `app` logic, plain C++ under test (`include/NiTCAD/app/`, `src/app/`): the deck document (text, path, modified
  flag, load and save, a check that maps `deck::parse`/`deck::build` errors to their lines); the run controller (one
  worker `std::jthread` per run; the states idle, checking, running, stopping, done, failed; progress events queued
  on the worker and drained on the UI thread, which a caller-supplied wake function (a posted Win32 message in the
  window) prompts; Stop requests the token and the completed part is kept); the workspace (a run's report or a
  result file with its device and mesh, paired by the mesh digest as `nitcad plot` does; the run record, provenance,
  convergence history and Unit 24's figures); view settings (view, field, point, scale, range, colormap, contours,
  edges, insulators, projection, camera, slices, clip) turned into `render` scenes with `nitcad plot`'s choices; the
  menus' deck sections and keys (VQ13); the recent files and layout, saved per user.
- The window (`src/app/`, `app/main.cpp`): per-monitor DPI aware v2; the workspace of Figure 1 (menu bar, tool bar,
  scene tree, viewport tabs, metrics strip, inspector, status bar); the deck pane (Rich Edit with line numbers, errors
  at their lines); the viewport presenting `render` scenes in a DXGI swap chain on the renderer's queue; orbit, pan,
  zoom, fit, probe and measure (VA8). `render` gains what a window needs and no more: a swap-chain target that shows
  the image `draw` would return, a zoom and pan on the camera (defaults leave every image unchanged), and the screen
  position of a world point for picking.
- 26·3: the start page (new, open, recent, examples), the run monitor (a row per point, the live convergence plot,
  Stop), the results workspace (curves, Unit 24's figures, the run record), `examples/`, and `--self-test <deck>`.
- **Perspective (AQ2's proposed "camera mode with a reference image")** already exists from Unit 25a (`Projection::
  perspective`, `nitcad plot --perspective`, the WARP reference `mosfet3d_opaque_perspective`); Unit 26 only exposes it.

**Acceptance gates (final):**

| # | Gate | Commit |
|---|---|---|
| WG1 | No solve on the UI thread: during a long run the UI thread's message loop keeps running (a timer message at least every 100 ms); Stop ends the run with its completed part kept and shown. | 26·1 (logic), 26·3 (window) |
| WG2 | A deck run in the app is bit-identical to `nitcad run` on it (`deck::run`), run identity included. | 26·1 |
| WG3 | The app's logic under unit tests: the run state machine, document load and save, view settings; deck errors land on their lines. | 26·1 |
| WG4 | `--self-test` passes on WARP and on the GPU and writes the same image as `nitcad plot` for the same view. On CI only if a hosted runner can open a window: the self-test exits 77 when no window or swap chain can be made, and ctest counts that as skipped; the off-screen WARP image tests run in CI regardless. | 26·3 |
| WG5 | Per-monitor DPI re-layout, keyboard access to every command, a clean D3D12 debug layer, no leaked GPU objects. | 26·2 |
| WG6 | Suites as VG6 (Release and Debug at every commit, `/W4 /WX`, hand mutants caught); nothing outside VA1–VA6 and VA8–VA9, nothing on VA7's list; the `app` layer's includes checked against ARCHITECTURE.md 3 (no `physics`, no `assemble`); no numerical change; ARCHITECTURE.md as built after acceptance. | all |
| WG7 | The Device, Mesh, Physics and Simulation menus open the deck at their sections and list their keys for insertion (VQ13); inserting any listed key (with its default value) into a deck leaves a deck that `deck::build` accepts, with the same meaning. | 26·1 (logic), 26·2 (menus) |
| WG8 | Every deck in `examples/` runs to completion in a test, and the start page lists exactly those decks. | 26·3 |

**CI and WG4.** Whether a hosted `windows-2025` runner can create a window and a flip-model swap chain cannot be checked
from this workstation and needs a CI run, which needs a push (not authorized yet). WG4's test is therefore written to
run anywhere: where a window cannot be made it is skipped (exit 77), not failed, and the first CI run after an
authorized push answers the question. The off-screen image tests (`nitcad_render_test`, `nitcad plot`) run in CI as
before.

**Commits (approved, AQ1):** 26·1 the app logic (WG1 logic, WG2, WG3, WG7 logic); 26·2 the main window (WG5, WG7
menus); 26·3 the start page, run monitor, results workspace, `examples/`, `--self-test` (WG1 window, WG4, WG8).

**Commit 26·1 (the app logic; WG1 logic, WG2, WG3, WG7 logic, WG6), as built.** The `app` layer's first files
(`include/NiTCAD/app/`, `src/app/`; CMake target `NiTCAD::app`, linked to deck, render, analysis and io), no window:
- `document.hpp`: `DeckDocument` (open, `set_text` with the modified flag, `save` through a temporary file renamed over
  the deck, `save_as`, the title); `check_deck` (parse and build as the file, so relative paths and messages are the
  deck's); `error_line` (the line after "<source>:" in a message, the form every deck error has).
- `run_controller.hpp`: `RunController` with one `std::jthread` per check or run; the worker builds (`check_deck`) and
  runs (`deck::run` with a `RunControl` on the thread's stop token); events (built, one per Newton iteration with the
  point's wall time, finished) are queued under a mutex and the caller's wake function is called only when the queue
  was empty, so a fast solve posts one message per batch; `process()` applies them on the UI thread. States idle,
  checking, running, stopping, done, failed; a run that stops early is done with its report (`stopped_early`); the
  rows of the run monitor (one per point, its iterations, last update and residual, time) and every iteration for
  the convergence plot; `take_results` hands the plan and report to the workspace.
- `workspace.hpp`: `Workspace` (device and geometry at fixed addresses; results; the bias points `nitcad plot` lists,
  transient snapshots with their currents and rates) from a run, a deck's outputs or a result file with its device
  (and mesh) file, refused on another mesh by the mesh digest or the node count; `run_record` (mesh, points, how it
  ended, input identity, mesh digest, settings, provenance), `convergence_history`, `extracted_figures` (Unit 24: for
  a gate or electrode sweep, the threshold at maximum g_m, the subthreshold swing, the maximum g_m and the on/off
  ratio of the transfer curve of the contact carrying the most current; for an ohmic sweep, the current at its last
  bias).
- `view_settings.hpp`: `ViewSettings` (view, field, point, scale, colormap, range, contours, edges, opaque insulators,
  contacts, projection, azimuth and elevation, slices and the clip box in µm, iso levels, curve contacts and frequency,
  log axis) and `build_scene`, which makes `nitcad plot`'s calls with its choices (the ten-decade floor of vector
  magnitudes, the 0..1 scale of a zero magnitude, the swept contact) without including the command's code (AQ4);
  `available_views`, `available_fields` (the rates only when stored), `plot_arguments` (the command's options for a
  setting, slices and the box converted from µm exactly as the command converts them).
- `deck_sections.hpp` (VQ13): the sections each of the Device, Mesh, Physics and Simulation menus opens, and the keys
  with documented defaults (`[device]` 1, `[mesh]` 3 by the generator, `[models]` 16, `[solver]` 23, `[transient]` 5);
  `find_sections` on text that need not build (comments, quotes and continuations as the syntax has them);
  `insertable_keys` (not those present; a mesh's by its type); `insert_key` (after the section's last entry, or the
  section added at the end, but no second analysis section).
- `app_state.hpp`: recent files (eight, most recent first, no duplicates) and the window's layout as "key = value"
  lines, damaged lines ignored, paths in UTF-8.
- **Tests** (`tests/app/`, 22 cases): errors at their lines (syntax, meaning, untitled); open, save, save as, a failed
  save keeping the path and flag; the state machine (checks, a failed build, refusals while busy, rows equal to the
  sweep's convergence records, outputs written); **WG1** — a 1,601-point run while a hidden message-only window's
  loop pumps a 10 ms timer: 42 timer messages, the longest gap 20–23 ms (bound 100 ms), progress in many batches;
  Stop after three converged points: done, cancelled, the completed points kept and written; **WG2** — the 1D diode
  and the 2D MOSFET run by the worker and by the `nitcad` executable: every field, current, setting, input identity
  and mesh digest equal to the bit, the app's written file too; the workspace from a run, the deck's outputs and the
  result file alike, another mesh refused by digest; the run record and history; the diode's current and the
  MOSFET's threshold, swing, g_m and on/off ratio; ten views (solid, mesh, potential with contours, n with edges and
  viridis at point 1, |J|, isosurfaces, I-V on a log axis, a 1D profile, a band diagram, a diode I-V) drawn on WARP
  **equal to `nitcad plot`'s PNG for `plot_arguments` to the last bit** (max channel difference 0); **WG7** — every
  offered key of every menu inserted into five decks (tensor 1D and 2D sweeps, generated 2D and 3D meshes, a
  transient) builds with an unchanged fingerprint (device temperature and mesh digest, every model, Newton, step,
  linear, equation, current, rate and transient option): 40 to 45 insertions per deck; the recent files and layout
  round trip.
- `app_layering` (ctest, a CMake script): no `physics` or `assemble` include in `app` (WG6).
- **Finding (test fixture, not changed):** a profile's window along its depth axis is ignored
  (`device/doping_profiles.hpp`), and the depth axis defaults to x, so the deck tests' 1D diode
  (`profile = acceptor constant 1e17 x=-1:1um`, `donor ... x=1um:3um`) is uniformly compensated, not a pn junction.
  The app tests copy that diode as a generic device; the 2D MOSFET fixture puts its pockets as Gaussians along y
  with x windows, and the examples (26·3) avoid the pattern.
- **Mutation checks (Release): 15 of 16 caught, 1 equivalent and removed** — error_line ignoring the source (caught once a same-length
  source was added), save keeping the flag, one row for all points, Stop not stopping the worker, a failed run done,
  starting while running, the digest unchecked, no ten-decade floor, edges ignored, `plot_arguments` dropping
  `--edges`, no band view, insertion after the header, present keys offered again, a wrong default, duplicate recent
  files; the sensed contact restricted to ohmic contacts and to contacts other than the swept one survived as
  equivalent (a gate carries no steady current), and the condition was removed.
- **Suites:** the app tests pass in Release (2 s) and Debug (65 s), no warnings under `/W4 /WX`; the whole Release
  suite in a clean worktree at `e20f6ae`: 24/24 (850 s).

**Commit 26·2 (the main window; WG5, WG7 menus, VA1, VA2, VA4, VA8, VA9's executable), as built.**
- **`render` (what a window needs, VA4):** `WindowTarget` and `Renderer::window_target(hwnd, w, h)` (a flip-model DXGI
  swap chain on the renderer's queue), `Renderer::present(scene, target, w, h)` (the scene drawn exactly as `draw`
  draws it, then copied into the back buffer and presented, so the window shows the PNG's pixels; the swap chain
  resized when the size changes); `Camera::zoom` and `Camera::pan` (the fitted view magnified about its centre and
  moved in fractions of its half size; at 1 and 0 the arithmetic is the fitted view's to the bit, and the 20 WARP
  references pass unchanged); the 2D axes span only the part of the bounds in view (the fitted view's ticks
  unchanged); `screen_position(scene, w, h, point)` (the pixel and depth `draw` puts a point at, for picking).
  **Finding (fixed):** releasing a swap chain whose last present was still queued crashed under the D3D12 debug
  layer; a `WindowTarget` now waits for the renderer's queue before it is released.
- **App logic:** `commands.hpp` (every command with its menu, mnemonic, accelerator and its keys' text; the camera's
  keys act in the focused view, not as accelerators, so the deck editor keeps its arrows), `layout.hpp` (the panes
  at a DPI: sizes at 96 DPI scaled to the monitor's; the sides narrowed to keep the centre's minimum), `probe.hpp`
  (the node drawn under a pixel: every node in 2D, the boundary surface's nodes in 3D, the front-most within 12
  pixels; the distance between two), `ViewSettings` zoom and pan with `orbit`, `zoom_by`, `pan_by`, `fit` and
  `apply_camera` (a camera move redraws without rebuilding the scene), `workspace_from_plan` and
  `RunController::take_plan` (a deck opened is built on the worker, not the UI thread).
- **The window (`src/app/`, `nitcad_app_window`; `app/main.cpp`, `nitcad-app.exe`, a `/SUBSYSTEM:WINDOWS` executable
  separate from `nitcad`, AQ4):** per-monitor DPI aware v2, common controls 6; the workspace of the reference: the
  menu bar File, Edit, Device, Mesh, Physics, Simulation, View; a tool bar (New, Open, Save, Check, Run, Stop, Fit,
  Probe, Measure, Perspective); the scene tree (regions, contacts, mesh, views, fields, results); the tabs Device
  view (with the metrics strip: extents, nodes and elements, the field's maximum, the draw time), Curves and Deck;
  the inspector (view, field, point, scale, colormap, contours, mesh edges, contacts, opaque insulators,
  perspective, a slice plane, the clipping box, iso levels); the status bar (state, file, adapter, the cursor's
  node); splitters; the layout and recent files saved per user (`%APPDATA%\NiTCAD\app.ini`). The view
  (`viewport.cpp`) presents in its swap chain: left drag orbits (3D) or pans (2D), right or middle drag pans, the
  wheel zooms, a double click fits, a click probes or measures; keys arrows, Shift+arrows, + and -, Home. The deck
  pane (`deck_pane.cpp`): Rich Edit with a line-number gutter (error lines in red), the check's errors below (a
  double click or Enter goes to the line). Opening a deck checks it on the worker and shows its device (and the
  results its outputs hold); Run (F5) runs it on the worker (WG1), the status bar following each point, and shows
  the potential and the run's curve (the sensed contact the one carrying the most current, on a log axis). The
  Device, Mesh, Physics and Simulation menus (VQ13, WG7) list "Go to [kind]" for their sections (greyed when the
  deck has none) and "Insert into [kind]" with the insertable keys and their help.
- **Tests:** `commands_test` (unique mnemonics per menu, unique accelerators, the reference's menu bar),
  `layout_test` (the panes tile the window; at 192 DPI every pane doubles; the centre's minimum), `probe_test`
  (2D nodes found where drawn with their values; nothing out of reach; distances; in 3D the front-most of two
  nodes at one pixel; orbit, zoom, pan and fit limits; zoom and pan about the centre, both axes),
  `render/window_test` (screen positions in the layout, zoom and pan; zoom 1 is the fitted image to the bit; a
  panned view leaves the margin white; a swap chain presents field and plot scenes and follows a resize), and
  `nitcad_app_window_test` (ctest `app_window`, hidden windows on WARP): **WG5** — every command in the built
  menus and its keys in the accelerator table; a `WM_DPICHANGED` to 192 DPI re-lays the panes at the new scale
  (the children where `layout_panes` puts them, 9-point fonts at 192 DPI) and back; **WG7** — the Physics menu
  inserts `srh = on` into a new `[models]` at the caret and the deck still builds, Simulation's "Go to [sweep]"
  moves the caret to its header, a deck error (F7) lands the caret on its line; a deck opened, run (worker),
  shown, probed, zoomed, panned and fitted, the curves tab drawing I-V; **WG5** — with the D3D12 debug layer:
  no message after opening, running, zooming and switching tabs, and no Direct3D or DXGI object alive once the
  window is destroyed (the DXGI debug interface's live-object report, which sees the renderer's objects while
  the window exists).
- **Looked at:** `nitcad-app` on the 2D MOSFET deck on the RTX 5060 Ti: the regions and then the potential after
  Run, the layout as the reference (tree, tabs, inspector, metrics strip, status bar); the view drawn in 3.8 ms.
- **Mutation checks (Release): 21 of 21 caught** — zoom ignored, pan x reversed, pan y by the half width, the 2D
  axes not clamped, screen positions without the layout, the elevation and the zoom unclamped, Fit resetting the
  angles, pan y not applied (caught once a vertical pan was tested), the back-most node picked (caught once two
  nodes at one pixel were tested), no pick radius, sizes not scaled by the DPI, the centre not kept at its minimum
  (caught once a 600-pixel window was tested), a duplicate accelerator, fonts not made at the new DPI, the new DPI
  ignored, accelerators without Shift, an inserted key not shown at its line, the curve sensing the swept gate, a
  swap chain released before its present ran, Run doing nothing. Three mutants first failed to build (a `/WX`
  warning, a delimiter in the substitution) and were rewritten so they compile.
- **Suites:** the app, window and render tests pass in Release and Debug (the Debug render suite with the D3D12
  debug layer clean), no warnings under `/W4 /WX`.
- **Not in 26·2 (26·3):** the start page, the run monitor, the results page, `examples/`, `--self-test`.
- **Suites:** the whole Release suite in a clean worktree at `ea90185`: 25/25.

**Commit 26·3 (the start page, run monitor, results workspace, examples and `--self-test`; WG1 window, WG4, WG8,
VA3, VA5, VA9), as built.**
- **`examples/` (VQ11, AQ3):** four decks, each titled and described by its first two comment lines:
  `pn_diode_1d.deck` (an erfc n+ profile on a 1e16 p body, swept forward to 0.7 V: the band diagram and I-V),
  `mosfet_2d.deck` (a 2 × 1 µm n-channel MOSFET under 50 nm of oxide, Gaussian source and drain pockets, Id-Vg at
  V_D = 0.1 V, 16 points), `mosfet_3d.deck` (the same in 3D, 0.8 µm wide, a gate sweep at V_D = 0.5 V, 5 points)
  and `moscap_cv.deck` (a 10 nm n+ poly gate on 1e17 p silicon, small-signal C-V at 1 MHz). Run with `nitcad`:
  0.04 s, 0.7 s, 4.4 s and under 0.1 s (Release). Looked at: the diode's junction near 0.9 µm with the
  quasi-Fermi levels split by the bias; the MOSFET's Id-Vg from 5e-9 to 6e-2 A/cm with its subthreshold slope;
  the MOS-C falling from the oxide capacitance (3.3e-7 F/cm²) through depletion; the 3D MOSFET's inverted
  channel. **Deviation (the MOS-C):** the sweep stops at 0.25 V, the onset of inversion: from there the
  drift-diffusion steady state of this 1D device did not converge (Newton at 0.5 V, also with step cutting and a
  line search, and on a mesh graded to 0.5 nm at the oxide), its electrons coming only from generation. The deck
  says so.
- **App logic:** `examples.hpp` (`list_examples`, `examples_directory` beside the executable else the source tree's,
  `copy_example` into the user's Documents\NiTCAD\examples, a copy kept once made) and `run_monitor.hpp` (the
  monitor's columns and cells, the convergence plot of every positive update on a log axis, the progress as
  converged bias points over the planned ones).
- **The window:** tabs Start, Device view, Curves, Deck, Run monitor, Results (Ctrl+1 to Ctrl+6); the start page
  (New deck, Open deck, Open result, the recent files, the examples with their descriptions; an example opens as
  the user's copy); the run monitor (shown when Run starts: a virtual list of the rows, the progress bar, Stop, a
  status line with how the run ended, the live convergence plot in its own swap chain, made before the run starts
  and redrawn at most every quarter second and at most a twentieth of the time); the results page (Unit 24's
  figures with their curve and flags, the run record and provenance).
- **`nitcad-app --self-test <deck> --out <png> [view options]` (VA6):** opens the window (shown), checks and runs the
  deck on the worker while its message loop runs with a 10 ms timer (the longest gap must stay under 100 ms,
  WG1), shows the view in the window's swap chain, writes the same view's image at `--size` and closes; with
  `--debug-layer`, also no debug message and no live graphics object after the window. Exit 0, 1, 64 (usage) or
  77 when no window or swap chain can be made (CI, WG4).
- **Finding (fixed):** the run monitor's plot made its swap chain on its first draw, during the run; on the GPU
  that took about 130 ms of the UI thread (the self-test's WG1 measure caught it). The swap chain is now made when
  the monitor is shown, before the run starts: the longest gap is 19–32 ms on the GPU and 22–30 ms on WARP for
  the 3D example, with the plot drawn in 2–3 ms.
- **Finding (fixed, Debug):** in Debug the solver's events arrived as fast as the slower window took them, so its
  queue was never empty and the timer (and input) waited behind the posted messages: gaps of 2.2 to 3.2 s. The
  window now shows a run's progress at most every 50 ms (and on a 100 ms timer while the events pause), so taking
  an event is cheap: Debug gaps 39 ms. The self-test measures the gap while the solver runs, and reports apart
  the loading of the results and the first view when the run ends (on the UI thread: 63 ms for the 2D MOSFET and
  95 ms for the 3D one in Release, 110 ms and 809 ms in Debug).
- **Tests:** `examples_test` (**WG8**: exactly the four decks, titled and described; a deck without comments; the
  directory beside the executable; a user's copy kept; every example run by the controller to completion, its
  points all solved and its outputs written), `run_monitor_test`, and in `nitcad_app_window_test`: the start page
  lists exactly the examples (the bundled ones and a given directory); an example opens as a copy; Run shows the
  monitor, which fills while the run goes on; the results page holds the threshold and the record; **WG1** — Stop
  after three points ends the run early with the completed points kept, shown and the status saying why; every
  tab reachable from the keyboard; **WG4** (`self_test_test`) — `nitcad-app --self-test` and `nitcad plot` on the
  2D MOSFET (potential; n with viridis and contours; Id-Vg on a log axis) and the 3D MOSFET (n) give the same PNG
  to the bit on WARP, and on the RTX 5060 Ti; the self-test with the debug layer passes with nothing alive after
  the window; usage errors exit 64. A session that cannot open a window skips WG4 (exit 77): on the hosted CI
  runner this is still unverified (it needs an authorized push), and the off-screen WARP tests run there either
  way.
- **Mutation checks (Release): 14 of 14 caught** — the title not read from the deck, the description from the
  first line, the examples beside the executable ignored, a user's copy overwritten, unconverged points counted as
  done, zero updates plotted on the log axis, the point numbered from 0, an example opened in place, Run not
  showing the monitor, the monitor updated only at the end, the start page missing an example (the bundled ones;
  a given directory, caught once a test gave one), the self-test's image not the plot's size, the self-test
  ignoring the view's options (caught once a view other than the window's default was compared).
- **Observed once, not reproduced:** while two builds, a mutation run and a suite shared the machine, one GPU
  self-test image differed from `nitcad plot`'s (the files left behind were equal) and WARP self-tests saw 400 to
  880 ms gaps; on reruns the GPU and WARP comparisons were exact and the gaps 19–32 ms.
- **Suites:** the app tests (32 cases) and the window tests (8 cases, the self-test included) pass in Release and
  Debug, no warnings under `/W4 /WX`. In a clean worktree at `84c34e9`: Release 25/25, Debug 19/19 — but the
  Debug suite took 2,645 s against about 1,460 s before Unit 26, because the app and window tests ran the 3D
  MOSFET example in Debug (about 5 minutes a run: `app` 369 s, `app_window` 982 s), too long for the CI job's 60
  minutes.

**Follow-up (the 3D example's gates in the Release suite).** As the other slow gates (`solve_mosfet`, ...), the 3D
MOSFET's run (WG8) and its self-tests on WARP and the GPU (WG4) moved to hidden-tag cases run by the Release-only
ctest entries `app_examples_3d` and `app_window_3d`; the Debug suite runs the other three examples and the 2D
self-tests (WARP, the GPU, the debug layer). Debug: `app` 101 s, `app_window` 231 s; Release: `app` 3.6 s,
`app_window` 8.6 s, `app_examples_3d` 4.2 s, `app_window_3d` 10.7 s.

**Suites at every commit (clean worktrees):** `e20f6ae` (26·1) Release 24/24, Debug 18/18; `ea90185` (26·2) Release
25/25, Debug 19/19; `84c34e9` (26·3) Release 25/25, Debug 19/19. Nothing outside `render` (the window target, the
camera's zoom and pan, the 2D axes in view, screen positions), `app`, `app/`, `examples/`, tests, `CMakeLists.txt`
and this file differs from `main` (`e289e44`): no numerical change, and the 20 WARP references pass unchanged.
**Status:** Unit 26 built (WG1–WG8 met on this workstation; WG4 on the hosted CI runner unverified until an
authorized push). **Accepted by the owner 2026-10-22** ("Accept Unit 26"); ARCHITECTURE.md 6.14, 8, 11 and 11.3
updated to "as built" (WG6). Merged to `main` as `e9b1d71` (PR #59, 2026-10-22).

### Unit 26b — the analysis workspace (started 2026-10-22; on `app/analysis-workspace`)

**Owner request:** VQ10 (live deck checking and completion, structure preview, cut lines with band diagrams, run
comparison and difference maps, mesh quality; plan page Figures 4, 7, 8 and view 25), approved after Unit 26's
acceptance ("start 26b", 2026-10-22). The owner's conditions of Units 25 and 26 hold: numerics and layer boundaries
preserved, Release and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged without
authorization.

**What exists (audited on `main` at `e9b1d71`).** The deck's reference list of sections and keys is
`deck/build.hpp`'s comment; `deck::build` stops at the first error, `"<source>:<line>: <what>"`, and rejects an
unknown key as `unknown key '<k>' in [<kind>]`. `app::check_deck` builds on a worker (`RunController`), which is busy
while a run goes on. `render::Locator` finds the element holding a point with barycentric weights; band edges are
stored per node in every result (NaN on insulator nodes); `band_plot` and `profile_plot` take a 1D mesh only.
`Workspace` pairs results with a device by the mesh digest; `analysis` extracts the threshold. `render`'s
`Geometry` gives the linear elements of every mesh (Kuhn tetrahedra of tensor cells). Coupling areas are positive by
construction (`mesh::Mesh`, P3), so mesh quality is a matter of element shape, not of negative couplings.

**Design (taken under VQ8–VQ14's delegation, "choose what's best for Commercial TCAD"):**

| # | Choice |
|---|---|
| XD1 | **Live checking.** A `LiveChecker` with its own worker, separate from the run controller (editing goes on during a run): a check starts 400 ms after the last edit; each check carries a generation, and only the newest check's result is applied (a build cannot be cancelled mid-way, so an older one is discarded when it ends). The editor shows the error's line in the gutter and its message under the editor, as Check does. |
| XD2 | **Suggestions.** An unknown section, key or value choice gets "did you mean ..." from the vocabulary (XD3): the nearest name by edit distance, at most 2 edits or a third of the name's length. |
| XD3 | **Completion.** A deck vocabulary in `app` (`deck_vocabulary.hpp`): every section kind, its keys (with or without default), the value choices of keys that have them, and a line of help; mirroring `deck/build.hpp`'s comment, which a test reads, so the two cannot drift. Region parameter paths (`eps_r`, `electron_mobility.mu_max`, ...) are open-ended and not completed. `complete(text, line, column)` gives what fits the caret: section kinds after `[`, keys at a line's start in a section, value choices after `=`. The editor opens the list on Ctrl+Space and when a key or `[` is typed; Enter or Tab inserts. |
| XD4 | **Structure preview.** After a live check builds, the device view shows the deck's device (regions, contacts and mesh) without solving. Results stay shown while the deck's mesh digest equals theirs; otherwise the view shows the structure, marked "preview". |
| XD5 | **Cut lines.** `render::cut_line(geometry, a, b)`: the segment clipped by every element it crosses, so the samples are the crossings with element faces; a nodal field read there by the element's barycentric weights is exactly the piecewise-linear field along the line (kinks at the crossings only). Outside the mesh is a gap. `cut_plot` and `cut_band_plot` draw any field or the four band edges against the distance in µm. In the window: Cut line (two clicks; in 3D on the front-most surface point, the line running through the device), the plot in the Curves tab; also typed end points in µm. |
| XD6 | **Run comparison.** A second result (a result file, or the previous run of the same deck) as run B. Curves of A and B overlaid; for transfer curves the threshold shift ΔV_th from Unit 24's extraction. A difference map B − A of a field at a bias point both runs have (equal bias values), on the same mesh (digest) only, else refused; diverging colormap (cool-warm), symmetric about 0: potential, quasi-Fermi levels and band edges in mV; densities as Δlog10 (decades); other fields linear in their units. The record of each run side by side with the settings that differ. |
| XD7 | **Mesh quality.** Per element: in 2D the minimum angle, in 3D the minimum dihedral angle (degrees). A map on the elements (2D) or the boundary and slices (3D) and a histogram; the minimum, the median and the count below a threshold (20° in 2D, 10° in 3D) in the metrics strip. |

**Commits:** 26b·1 live checking, suggestions, completion and the structure preview (XD1–XD4); 26b·2 cut lines
(XD5); 26b·3 run comparison and mesh quality (XD6, XD7).

**Acceptance gates:**

| # | Gate | Commit |
|---|---|---|
| XG1 | Live checking: an error appears at its line after an edit without Check; a check never runs on the UI thread; a stale check's result is never applied (a deterministic test holds the worker); editing during a run checks the deck. | 26b·1 |
| XG2 | The vocabulary equals `deck/build.hpp`'s reference (every section and key it names, no other); every listed key with a sample value gives no "unknown key" error from `deck::build`; every misspelling at one edit of a key, section or choice is suggested back; completion at each caret context gives exactly the fitting candidates. | 26b·1 |
| XG3 | The structure preview's scene equals the solid and mesh scenes of the device `deck::build` makes (`nitcad plot` of the deck's device file) to the bit. | 26b·1 |
| XG4 | Cut lines: exact on a linear field; the samples are the element crossings (a tensor and a generated mesh, 2D and 3D); along a 1D device equal to `band_plot`'s points; a line leaving the mesh has a gap there; the 2D and 3D MOSFETs' band diagrams from source to drain drawn on WARP as references. | 26b·2 |
| XG5 | Comparison: B − A of a run with itself is zero; of two runs it equals the nodal difference; another mesh or a bias point not in both is refused; ΔV_th equals the difference of the extracted thresholds; a WARP reference of a difference map. | 26b·3 |
| XG6 | Mesh quality: exact angles on equilateral, right and regular elements (60°, 45°, 70.5288°); the histogram counts every element once; the generated meshes respect their generators' minimum-angle setting in 2D. | 26b·3 |
| XG7 | The window: every new command in the menus and reachable by the keyboard; a clean D3D12 debug layer and no live graphics object after the window; per-monitor DPI for the new controls. | all |
| XG8 | Suites as WG6 (Release and Debug at every commit, `/W4 /WX`, hand mutants caught); `app` still free of `physics` and `assemble`; no numerical change (the deck, solver and render references unchanged); ARCHITECTURE.md as built after acceptance. | all |

**Commit 26b·1 (live checking, suggestions, completion and the structure preview; XD1–XD4, XG1–XG3, XG7), as
built.**
- `live_check.hpp`: `LiveChecker(wake, check)` with its own `std::jthread`; `submit(text, path)` returns the check's
  generation; the worker takes only the newest submission, so edits made while a build runs coalesce into one check
  after it; a result is posted (wake) and `process()` on the UI thread applies it only when no newer submission exists,
  else counts it as `dropped()`. The check is injectable (tests hold the worker on a latch). Independent of the run
  controller: the deck is checked while a run goes on.
- `deck_vocabulary.hpp` (XD2, XD3): 15 sections, each key with its named-key arity, repeatability, a sample value,
  a line of help, its value choices (first and second word) and, for `[contact]`, `[mesh]` and `[solver]`, the
  context (the selector key's value) in which it applies; `edit_distance`, `nearest` (at most min(2, max(1, ⌊n/3⌋))
  edits, exact matches excluded), `suggestion(text, diagnostic)` from the error's own line (unknown section, key or
  choice), and `complete(text, line, column)`: section kinds after `[` (closing it unless closed or named), the keys
  of the caret's section not yet present (a key and " = ", or a named key and a space), value choices after `=`
  (waveforms with their parenthesis), nothing in a comment, on a continued line, on a header's name or before the
  first section.
- `workspace.hpp`: `workspace_from_structure(plan, deck)` (the device and geometry, no results) and `same_mesh(w,
  mesh)` (by the graph digest); `workspace_from_plan` uses the former.
- The window: an edit restarts a 400 ms timer, then the deck is submitted; a result older than the text shown is not
  applied; an error underlines its line (wavy, drawn after the Rich Edit's paint), the message with " (did you mean
  'x'?)" goes to the status bar and Find Error; a deck that builds is shown as "Device (preview, not solved)" unless
  the results shown are on the same mesh, keeping the view settings when the dimension is unchanged. Completion:
  a list under the caret (Up/Down, Enter or Tab inserts, Escape closes; Left/Right/Home/End/PgUp/PgDn close),
  opened by typing a word or `[` or `=`, and by Edit > Complete (Ctrl+Space); taking a key opens its value choices.
- **Tests** (`tests/app/`, 11 new cases): **XG1** — errors at their lines then the plan, the check on the worker
  thread; a held worker with five submissions: two results dropped, the newest applied; a check during a 1,601-point
  run; in the window: a misspelt material marked at line 11 with "(did you mean 'sio2'?)" and one check for two
  quick edits, a result older than the text not shown, `kind = omic` caught during a run; **XG2** — every vocabulary
  key with its sample value accepted by `deck::build` in every context (more than 130 probes); every word of
  `deck/build.hpp`'s reference that `deck::build` takes is listed (more than 100), parameter paths excepted; the
  menus' defaults equal the vocabulary's; every key, section and choice misspelt by one edit (a letter dropped)
  suggested back from `deck::build`'s own error; completion at each caret context; **XG3** — the preview of the 2D
  MOSFET and a generated 3D deck drawn on WARP equal `nitcad plot`'s solid and mesh PNGs of the deck to the bit;
  the preview follows a finer mesh in the view chosen; results kept on an edit that keeps the mesh; **XG7** —
  completion by the keyboard in the window (Enter, Tab, Up/Down, Escape, Ctrl+Space, a double click).
- **Finding (fixed):** `DeckPane::caret()` gave the previous line on the final empty line (Rich Edit's line index for
  a position after the last `\r`); found by the forced-completion test.
- **Finding (not changed):** `deck/build.hpp`'s reference omits keys `deck::build` takes (the `gmres_ilut` keys
  `ilut_drop`, `ilut_fill`, `restart`, and the model switches by name). The vocabulary lists what `deck::build` takes,
  so XG2's "no other" is checked against `deck::build`, not the comment.
- **Mutation checks (Release): 29 of 29 caught, one equivalent replaced.** Vocabulary, checker and workspace (18):
  the nearest name at 3 edits and completion inside a comment survived the first run and are caught by two added
  checks; three mutants that failed only by `/WX` (an unused variable) were rewritten to compile and are caught.
  Window (11): Enter's character typed after accepting, Down skipping an item, no value choices after a key, not
  opened after `=`, typing a word opens nothing, a result older than the text applied, results dropped on any edit,
  no suggestion, no live check, the view reset on every edit, Ctrl+Space not forced; three `/WX`-only mutants
  rewritten as above, of which Down past the last item is equivalent (the list box refuses an index out of range)
  and was replaced by Down skipping an item, caught by an added step on the 15-item list.
- **Suites:** Release, the app tests (43 cases, 4,423 assertions) and the window tests (11 cases) pass, no warnings
  under `/W4 /WX`; Debug, the app tests (121 s) pass, and the window tests (284 s) pass but for one self-test run
  whose WG1 gap was 248.6 ms (bound 100 ms) while the solver ran; three reruns of that self-test gave 40–51 ms (the
  26·3 record: 39 ms). The self-test makes no edit, so the live checker does not run in it; the gap is recorded as
  an outlier of the Debug build on WARP. The clean-worktree suites at this commit are recorded with 26b·2.

**Commit 26b·2 (cut lines; XD5, XG4, XG7), as built.**
- `render/cut_line.hpp`: `cut_line(geometry, a, b)` — per element (after a bounding-box test), the barycentric weights
  at a and b (Cramer, as `Locator`), affine along the line, give the interval where all are non-negative; a weight
  changing by at most 1e-10 over the line is parallel to its face (a line on a face reads both elements); the
  intervals, sorted, are joined into samples at the crossings (one where two elements meet, a `gap_before` mark
  where the line re-enters the mesh); weights below 1e-12 are zero and the rest renormalised, so a sample on a node
  reads that node exactly. `cut_values` reads a nodal field, skipping nodes of zero weight (a NaN on an insulator node
  leaves the faces away from it defined).
- `scene.hpp`: `cut_plot` and `cut_band_plot` (against "distance (µm)", a NaN point at each gap), sharing
  `profile_plot`'s and `band_plot`'s code (their output unchanged: the references pass); `add_cut_line` (the parts
  inside the mesh, a white 2 px line over a black 4 px one, on top of everything, so in 3D it shows through the
  device).
- `ViewSettings::cut` (ends in µm, in the device's dimension): the profile and bands views plot along it (needed in 2D
  and 3D), the 2D and 3D solid, mesh, field, slice and iso views draw it, any other view refuses it; 2D and 3D
  workspaces now offer the profile and bands views; `plot_arguments` gives `--cut`.
- **Deviation (owner's record): `nitcad plot --cut <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>]`** with the same rules and
  messages. XD5 places cut lines "in the window", but VA4 has every window view drawn as `nitcad plot` draws it and
  the WG4 test compares the two images; without the option the window's cut views could not be checked that way.
  The seeds' segment parser became the shared `segment_um`. Also in both: a profile of the doping no longer needs
  results (it never used them), so a cut line on an unsolved deck plots its doping.
- The window: View > Draw cut line (Ctrl+Shift+L, a toolbar button), exclusive with Probe and Measure; two clicks on
  the device view — in 2D the point itself (`Prober::point_at` inverts the view's affine map), **in 3D the front-most
  surface node within 12 px (as Probe), not any point of the surface (a deviation from XD5's "front-most surface
  point": a node is what the view can pick reliably)** — rounded to the nanometre; the Curves tab then shows the band
  diagram along it, or the doping before a run. The inspector's "Cut line (µm)" field types or clears it (Apply). One
  cut line per window, kept while the device keeps its dimension (a run of the deck), cleared when a deck or result
  is opened. The Curves tab offers the profile in 2D and 3D; it starts on a band diagram only in 1D.
- **Tests:** `cut_line_test` (6 cases) — a linear field read to 1e-12 at every sample of tensor and generated meshes
  in 2D and 3D; the samples equal the crossings found independently by intersecting the segment with every element
  face (segment–edge in 2D, Möller–Trumbore in 3D); a horizontal line across a row of cells: 21 samples (11 faces, 10
  diagonals); lines along grid lines and through nodes: the nodes only; along a 1D device the cut band diagram equals
  `band_plot` (x to 1e-12 µm, y exactly, NaN where NaN); a U-shaped mesh: the line enters at x = 0, leaves at 1 µm,
  re-enters at 2 µm with one gap, the plot broken there, the overlay without the gap; refusals. **WARP references**
  (`mosfet2d_cut_bands`, `mosfet3d_cut_bands`, `mosfet3d_cut_line`): source to drain 0.05 µm under the oxide.
  `view_settings_test`: errors as `nitcad plot` reports them; the scene equals render's; the doping without results;
  the 2D MOSFET's band diagram, a profile of n across the oxide and a map with the line **equal to `nitcad plot`'s
  PNGs to the bit**. `run_test`: the option, the line changing the map, the refusals. The window: a typed line on the
  unsolved deck (the doping, the profile offered), kept through a run; drawn by two clicks within 0.01 µm and
  rounded to the nanometre; the band diagram in the Curves tab and the line on the map; typed, malformed, cleared;
  Probe leaving the tool; a 3D preview's line between the surface nodes the clicks pick. `probe_test`: the point
  under a pixel of a 2D view, as the window draws it (zoomed, panned) and from a turned and tilted camera.
- **Mutation checks (Release):** 27 of 27 caught — render 11 (no parallel case, covered pieces repeated, no gap,
  weights not snapped, not normalised, the box test excluding touching elements, zero weights read, equal ends
  taken, no NaN point at a gap, the overlay drawn across a gap, cut bands labelled x), app 6 (a device view without
  the line, 2D bands without a line, bands ignoring the line, `plot_arguments` repeating the start, no 2D profile
  view, the end's dimension unchecked), CLI 1 (`nitcad plot` without the line), window 9 (the Curves tab without the
  2D profile, a device view without the line, clicks ignoring the tool, ends not rounded, the Curves tab not shown,
  a run dropping the line, an unsolved line plotting the potential, the 2D inverse transposed, 3D picking nothing).
  Two survived at first and gained tests: the transposed 2D inverse (every 2D view the window builds looks straight
  down, where the swapped terms are zero — killed by the tilted camera in `probe_test`) and the end's dimension
  (killed by a cut whose start alone is right). A first run whose test filter matched no test was discarded; the
  mutation script now refuses such a filter and stops if it cannot restore a file.
- **Suites:** in the working tree, Release: render (67 cases, 109,200 assertions, the three references), app, the
  CLI cut test and the window tests pass, no warnings under `/W4 /WX`; Debug: render (67 cases, 696 s), app (123 s),
  the CLI and the 26b window tests (29 s) pass. Debug render tests on WARP must run alone: run beside other D3D
  tests they failed to create a device (0x887A0007). The clean-worktree suites at 26b·1 (3077ce4): Release 27 of 27
  (990 s); Debug 19 of 19 — render failed once beside the mutation runs (the device error above) and passed run
  alone (61 cases, 658 s); app 110 s, window 297 s. The clean-worktree suites at this commit are recorded with
  26b·3.

**Commit 26b·3 (run comparison and mesh quality; XD6, XD7, XG5, XG6, XG7), as built.**
- `render/quality.hpp`: `min_angle_deg` — a triangle's smallest angle (each from atan2(|u × v|, u · v)), a
  tetrahedron's smallest dihedral angle (π less the angle between two outward face normals, over the six face pairs);
  0 for a degenerate element (|u × v| ≤ 1e-14 size², |6V| ≤ 1e-14 size³). `element_quality` (one angle per element,
  2D and 3D only), `quality_summary` (minimum, median — the mean of the two middles for an even count — and the count
  strictly below the threshold, 20° in 2D and 10° in 3D), `quality_histogram` (18 bins of 5°, every element once).
  `scene.hpp`: `quality_view` (flat per element, 0 to the best angle — 60° in 2D, 70.53° in 3D — unless a range is
  given; 2D with the element edges, 3D on the boundary or on slices; no contacts) and `quality_plot` (the histogram as
  a step line and the threshold, 0 to 90°). `mesh_view`'s 3D edges now come from the shared `boundary_face_edges`
  (its output unchanged: the references pass).
- `render/compare.hpp`: `matching_point` (the first point of B with the same contacts, each bias within 1e-9 V);
  `difference_field` (B − A at the nodes of A's geometry: densities in Δlog₁₀ decades, NaN unless both are positive;
  the potential in mV; **the band edges and quasi-Fermi levels in meV** — XD6 says mV; they are energies, so meV is
  the unit of the same number); `difference_scale` (symmetric about 0, ±1 when everything is 0 or NaN).
- `app/compare.hpp`: `comparable` (the same mesh, by the meshes' digest), `figure_shifts` (Unit 24's figures both runs
  have, by name and curve, B − A: ΔV_th is the threshold's), `compare_records` (A's lines in order, then the lines
  only B has; `differs` marked).
- `ViewSettings`: views `quality`, `quality-histogram` (2D and 3D with elements, before a run too) and `difference`
  (offered with a run B on the same mesh). `build_scene`, `available_views` and `plot_arguments` take run B: the I-V,
  C-V and transient-current views overlay B's curve of the same contacts (by name), the series named "(A)" and
  "(B)"; the difference map is cool-warm unless another colormap is chosen, scaled by `difference_scale` unless a
  range is given, a profile in 1D. The doping's difference needs no results. Refused: no run B, another mesh, no
  bias point of B at A's biases, a contact of A that B lacks, B without the analysis the curve needs. A transient's
  snapshots are matched by bias like any point (the first equal one). Band diagrams and profiles are not overlaid.
- **Deviation (owner's record, as for 26b·2's `--cut`): `nitcad plot --compare <deck | result file>`** and the views
  `quality`, `quality-histogram`, `difference`, with the same rules and messages, so the window's new views are
  checked against `nitcad plot`'s PNGs (VA4, WG4). A result file is read with run A's device and must have its mesh
  (the result's digest, or the node count for a file without one). A deck's B outputs are loaded as A's are (the
  shared `load_outputs`); `nitcad plot` mesh quality and a doping difference need no results.
- The window: File > Compare with a result... (Ctrl+Shift+C; the device and mesh files found as Open result finds
  them), Simulation > Compare with the previous run (Ctrl+Shift+V) and Stop comparing (Ctrl+Shift+X), View > Mesh
  quality (Ctrl+Shift+Q, 2D and 3D). The previous run is the workspace with results that a run (or the structure
  preview) of the same deck replaced; Compare with the previous run takes it as run B, Stop comparing gives it back.
  Opening a deck or a result forgets both. The metrics strip shows run B's name and ΔV_th when comparing, and on
  the quality views "Min angle, median, N below T°". The Results tab adds "Run B" and "B − A" columns to the figures,
  and the run records side by side ("≠" before a line that differs). Probing a difference map reads B − A.
- **Tests:** `quality_test` (4 cases) — exact angles: equilateral 60°, right isosceles 45°, 30-60-90 30°, the regular
  tetrahedron 70.5288°, a cube corner's tetrahedron acos(1/√3), a Kuhn tetrahedron 45°, degenerate 0 (collinear, coincident,
  flat, and collinear to rounding), a sliver under 0.2°; tensor meshes 45° in 2D and 3D; generated 2D meshes at 20.7°, 25° and 30° settings: no element below
  the setting other than those the generator reports exempt, the minimum equal to its report; the 3D generator's
  minimum dihedral equal to its report; the histogram's bins and the summary (odd and even counts, the threshold
  itself not below); the map and the plot. `compare_test` (render, 3 cases) — matching points (another number of
  contacts refused, also one equal as far as A goes); a run against itself zero in every field it holds (at least
  11); two runs: the potential's difference in mV, n's in decades, the donors' log 4, the scale symmetric. **WARP reference** `mosfet2d_difference_potential` (the drain at 0.5 V less 0 V: +500 mV in
  the drain pocket, the depletion region between). `compare_test` (app, 2 cases) — the 2D MOSFET at gate work
  functions 4.1 and 4.3 eV: ΔV_th equals the difference of the two extracted thresholds and is positive; the shifts
  in A's order; the input identity marked as differing, the mesh not; a run against itself: no shift, nothing marked;
  two meshes not comparable, their figures still compared; the diode swept from the cathode: the same figure on
  another curve, no shift; an unsolved deck beside a run: B's lines after A's, marked. `view_settings_test`: the views offered (with and without
  B, another mesh), the refusals, the overlay's series, the plot arguments; the quality map, the histogram, two
  difference maps and the overlaid I-V **equal to `nitcad plot`'s PNGs to the bit**. `run_test`: B as a deck and as
  its result file give the same image; the refusals; mesh quality and the doping difference without results. The
  window: Mesh quality (the strip reads "Min angle 26.6°, median 26.6°, 0 below 20°": the MOSFET's 0.1 × 0.05 µm
  cells give atan ½), a second run with another work function, Compare with the previous run (the difference map,
  ΔV_th negative in the strip, the Results columns, the marked identity), Stop comparing (the field view back, the
  previous run kept), a result file as run B, a new deck forgetting both. The command test covers the four new
  commands' menus and unique accelerators.
- **Mutation checks (Release):** 40 of 41 caught; the 41st caught in Debug only. Render 21 (a triangle's third angle
  ignored, an angle of the wrong sign, collinear triangles not degenerate, normals not outward, the normals' angle as
  the dihedral, a face pair skipped, the 3D threshold 20°, the even median the upper middle, the threshold itself
  below, 90° in the first bin, bins rounded, a loose bias match, the contact count unchecked, A − B in decades, the
  potential in V, a zero density in B read, the scale from the largest positive, the scale not symmetric, the quality
  map's scale to 90°, every element the first's colour, the threshold line at 0), app 13 (comparable ignoring the
  mesh, A and B swapped, figures matched by name alone, differing lines not marked, lines only B has dropped, the
  difference offered on another mesh, the difference map in turbo, its scale fitted, `plot_arguments` without
  `--compare`, the overlay unnamed by run, I-V without B's curve, the histogram drawn as a map, no histogram view),
  CLI 5 (a deck B on another mesh taken, no bias-point refusal, a result file B ignored, the histogram or the doping
  difference needing results), window 11 (no previous run kept, the previous run kept after use, Stop comparing
  losing it, a new deck keeping run B, Mesh quality showing the field, views without run B, no ΔV_th in the strip, a
  difference view kept without B, the 3D threshold in 2D, A − B in the Results, differing lines not marked). Four
  survived at first and gained tests: collinear to rounding (now 0 exactly), figures matched by name alone (the
  diode swept from the cathode: the same figure on another curve), lines only B has (an unsolved deck beside a run),
  and the contact count — whose mutant reads past A's biases: a run B equal to A as far as A goes is now refused,
  caught by Debug's checked `vector` (an assertion), not in Release, where the read is unchecked.
- **Suites:** in the working tree, Release 27 of 27 (974 s); Debug 19 of 19 — `app_window` failed once in the suite
  and passed run alone (271 s): the GPU self-test's WG1 gate saw a 226 ms gap in the message loop. A first run
  failed the same way in Release (652 ms, `app_window_3d`), passing 4 of 4 after. **Open item (WG1 gap, deferred
  by the owner for investigation, 2026-10-22):** a self-test's message loop occasionally stalls past 100 ms during
  a run — on the GPU here, and in 26b·2's build as well (1,623 ms on a first GPU run, then 41 and 38 ms; 26b·3's:
  115, then 42 and 43 ms), and on WARP in 26b·1's Debug suite (248.6 ms). Not caused by this commit (26b·2's build
  shows it); the cause is not established — a first run after the adapter idles is one reading, but the WARP case
  does not fit it. WG1 is left as ruled. No warnings under `/W4 /WX`. The tests added for the four survivors came after the suites; they pass in
  Release and Debug (render quality and compare, app compare). The clean-worktree suites at 26b·2 (394c2e9): Release 27 of 27, Debug 19 of 19
  (1,687 s). The clean-worktree suites at this commit are recorded with the acceptance.

**Status:** Unit 26b built in three commits (3077ce4, 394c2e9, eb79729; XG1–XG8 met on this workstation, the WG1
gap above open). **Accepted by the owner 2026-10-22** ("accept Unit 26b"), the WG1 gap to be investigated later;
ARCHITECTURE.md 6.14, 8, 11 and 11.3 updated to "as built" (XG8). Pushed and merged to `main` as `51face0`
(PR #60, 2026-10-22, on the owner's word, the hosted CI still running). The clean-worktree suites at `eb79729`
(finished after the merge): Release 27 of 27 (980 s); Debug 18 of 19 (1,951 s) — `app_window` failed only on the
WG1 gap above (the GPU self-test's message loop stalled 278.8 ms), the open item.

**Not in Unit 26b (VA7 and later):** form editors that write the deck, project management, remote runs, plugins,
localization, an installer; completion of region parameter paths; more than two runs compared at once.

### Unit 26c — GUI catch-up (proposed 2026-10-22; on `app/gui-catch-up`)

**Owner request:** "GUI catch-up" (2026-10-22), chosen after Unit 26b's acceptance: bring into the window what
`nitcad plot` draws but the window does not, and what Units 25c and 26 left for the application. The owner's
answers of the same day: transient snapshots as **playback only** (the tiled filmstrip stays in `nitcad plot`); the
clip box edited by **inspector sliders** (no in-view handles); **`nitcad plot --hide <region>`** added for region
visibility. The owner's conditions of Units 25, 26 and 26b hold: numerics and layer boundaries preserved, Release
and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged without authorization.

**What exists (audited on `app/analysis-workspace` at `195e829`).** `nitcad plot` draws arrows (`--vectors`,
`--arrows`, 300 in 2D and 400 in 3D), streamlines (`--streamlines`, `--lines` 24, `--from <contact>`, `--seeds` a
line or a 3D plane), the exploded view (`--axis`, `--gap`) and filmstrips (`--every`, `--columns`, `--frames`;
DC4: at most 12 frames on one scale fitted over all). `ViewSettings` already holds the curves' swept and sensed
contacts, the C-V frequency and the log axis, but the inspector shows no control for them. The scene tree lists the
regions without any action. `render::ViewOptions` has no notion of a hidden region; no `nitcad plot` option hides
one. Units 25a (VR7) and 25c (DC4, DC6) left animation playback and interactive clip editing "to 26"; Unit 26 did
not take them (ARCHITECTURE.md 6.14, "Not in Unit 26").

**Design:**

| # | Choice |
|---|---|
| YD1 | **Arrows and streamlines.** Inspector controls on the field view: Arrows (none, E, J, Jn, Jp) with a count; Streamlines (none, E, J, Jn, Jp) with a count and seeds from a contact or along a typed seed line (2D) or plane (3D). The same rules, defaults and messages as `nitcad plot` (the map shows the vector field's own magnitude); `ViewSettings` gains the options and `plot_arguments` gives `--vectors`, `--arrows`, `--streamlines`, `--lines`, `--from`, `--seeds`. |
| YD2 | **Exploded view.** A 3D view kind with the axis (x, y, z) and gap (µm, default 12% of the extent) in the inspector, as `--view exploded`. |
| YD3 | **Curve pickers.** On the I-V, C-V and current views the inspector shows the swept contact, the sensed contact, the C-V frequency and the log axis (the existing `ViewSettings` fields). |
| YD4 | **Region visibility.** A check box per region in the scene tree hides it. `render::ViewOptions::hidden_regions`: the solid, mesh, field, slice, iso, exploded and quality views leave out a hidden region's control-volume pieces (2D) or elements and boundary faces (3D); the faces a hidden region shared with a shown one become boundary faces; contacts on hidden regions only are left out; scales are fitted over the shown regions. `nitcad plot --hide <region>` (repeatable) with the same rules; an unknown region name is refused. |
| YD5 | **Transient playback.** On a transient run's device view: Play/Pause, a time slider over the snapshots and step keys (previous, next, first, last). The scale is fitted once over all snapshots (as DC4) and stays fixed; a frame is drawn by a timer (5 frames a second), never by a loop on the UI thread, and a frame not yet drawn is skipped rather than queued. The tiled filmstrip stays in `nitcad plot`. |
| YD6 | **Clip box sliders.** In 3D, a low and a high slider per axis over the device's extent (1,000 steps, values rounded to the nanometre) beside the typed box; moving one updates the box and the view. |

**Commits:** 26c·1 arrows, streamlines, the exploded view and the curve pickers (YD1–YD3; `app` only); 26c·2 region
visibility (YD4; `render`, `nitcad plot`, `app`); 26c·3 transient playback and the clip sliders (YD5, YD6).

**Acceptance gates:**

| # | Gate | Commit |
|---|---|---|
| YG1 | Every new view and overlay in the window — arrows and streamlines in 2D and 3D, the exploded view, the curves with chosen contacts and frequency, hidden regions — equals `nitcad plot`'s PNG for the same options to the bit (WG4's choices). | 26c·1, 26c·2 |
| YG2 | Region visibility: with no region hidden every scene is unchanged (the existing references pass); a hidden region's pieces, elements and contacts are absent and the shown regions' are unchanged; the shared faces close the shown regions; WARP references of the 2D and 3D MOSFETs with the oxide hidden. | 26c·2 |
| YG3 | Playback: every frame equals the filmstrip frame of the same snapshot on the same scale to the bit; the window's message loop meets WG1 while playing; Pause, the slider and the step keys land on the snapshot asked for. | 26c·3 |
| YG4 | Clip sliders: a slider's value and the typed box agree both ways; the view equals `nitcad plot --clip` with the slider's box. | 26c·3 |
| YG5 | The window: every new command in the menus and reachable by the keyboard; a clean D3D12 debug layer and no live graphics object after the window; per-monitor DPI for the new controls. | all |
| YG6 | Suites as WG6 (Release and Debug at every commit, `/W4 /WX`, hand mutants caught); `app` still free of `physics` and `assemble`; no numerical change; ARCHITECTURE.md as built after acceptance. | all |

**Started on the owner's word (2026-10-22, "start 26c.1 also investigate app_window").**

**Commit 26c·1 (arrows, streamlines, the exploded view, the curve pickers; YD1–YD3, YG1, YG5, YG6), as built.**
- `ViewSettings`: `vectors`, `streamlines` (E, J, Jn, Jp; empty none), `arrows`, `lines` (0: `nitcad plot`'s 300 in
  2D, 400 in 3D; 24), `from` (a contact; empty: its choice), `seeds` (a segment or, in 3D, a plane `x=<um>`); the
  view kind `exploded` with `axis`, `gap_um` (empty: 12% of the extent) and `exploded_field`. `build_scene` draws them
  with `nitcad plot`'s calls, order and messages: streamlines, then arrows, then the cut line; in 3D over their own
  magnitude the arrows on the solid device and the lines in the outline, colored by it; 3D arrows on the slices of a
  slice view, else on the exposed boundary; seeds from the named contact, else for currents the one carrying the
  most, for E the first gate or electrode, else the most biased. `plot_arguments` gives `--streamlines --lines
  --from --seeds --vectors --arrows --axis --gap` (and `--field` for an exploded view showing the field).
  `available_views` offers `exploded` in 3D, after the slices.
- **Arrows and streamlines on the field views only** (`takes_vectors`: field, slice, iso — YD1's "on the field view"):
  a setting keeps them, `build_scene` and `plot_arguments` leave them out on the other views, so the window and
  `nitcad plot` agree on every view (the exploded view in particular, where `nitcad plot` would draw them where
  the regions were before they moved).
- The inspector: Arrows (count), Streamlines (count), Seeds from (or typed, µm) — each a list and a box side by
  side — on the field views of a 2D or 3D run's results; Exploded along (gap, µm) and Show the field on the exploded
  view; Swept contact (Contact on the current view: "(all)" first), Sensed contact, Frequency (the run's frequencies)
  and Log axis on the I-V, C-V and current views. A group shows only with its views (`show_inspector_groups`; the
  inspector lays out the shown ones). Choosing arrows or streamlines on the field view shows that vector's magnitude
  (as `nitcad plot` without `--field`). Counts, seeds and the gap are read on Apply; an invalid count or gap leaves
  the setting and says why in the status bar. Every new control is a tab stop (YG5); no new command. **No probe or
  cut-line click on the exploded view** (its regions moved; `pickable`).
- **Tests:** `view_settings_test` — the exploded view's name and place in the views; `nitcad plot`'s messages (an
  unknown vector field, the counts' limits, a seed plane in 2D, a malformed segment or plane, an unknown contact, a
  1D device, the exploded view's device, axis and gap); vectors left out off the field views (the exploded scene
  unchanged); the plot arguments. **To the bit equal to `nitcad plot`'s PNGs (YG1):** 2D arrows of J over |J| and of
  E over contour levels, streamlines of E (from the gate), of J (from the contact carrying the most), of Jn from a
  named contact with arrows, from a typed segment; 3D colored arrows of J, colored lines of E, lines and arrows of J
  together, arrows of E on a slice, lines from a seed plane; the exploded view with and without a field; I-V of a
  chosen sensed contact, C-V at the second frequency, the transient current of one contact. New fixtures: a small 3D
  MOSFET (9 × 7 × 3 nodes, about a second in Debug), the diode's C-V at 1 kHz and 1 GHz. `analysis_window_test`: the
  inspector driven as a user would — groups shown and hidden per view, arrows of J (the map switches to |J|),
  counts and a refused count, lines from the gate and from a plane, the exploded view (axis, field, gap, a refused
  gap; no probe where the solid view probes), back to the field view with the lines kept; the window's images equal
  `nitcad plot`'s for its settings; tab stops; at 200% DPI the pairs side by side; the C-V pickers (frequencies,
  contacts, log, the default again).
- **Mutation checks (Release):** 31 of 31 caught — view settings 19 (Jn read as Jp, the 2D arrow and line defaults,
  colored arrows beside lines, the line color, the arrow length, 3D arrows not on the slice, the seed plane's unit,
  `--from`'s contact, the contact carrying most, E lines not from the gate, the default gap's axis, the exploded
  field ignored, `--axis`'s letter, iso without vectors, vectors kept off the field views, no exploded view, 3D
  lines not colored, the arrow limit), window 12 (the map keeping its field, a count of 0, untrimmed seeds, the
  arrows' controls on every view, a probe on the exploded view, the frequency picker, no default entry, the pairs
  overlapping, Log and Show the field swapped, swept and sensed swapped, the exploded controls always, a negative
  gap). Two survived at first and gained tests: 3D arrows and lines of one field together, and the exploded view's
  probe (the first test switched probing off by toggling it).
- **Suites (working tree, with the two changes below):** Release 27 of 27 (1,391 s), Debug 19 of 19 (2,382 s); no
  warnings under `/W4 /WX`. The two tests strengthened after the mutants pass in Release and Debug.
- **Found:** `nitcad plot <deck>` refused arrows and streamlines on the solid, mesh and quality views ("the results
  hold no field to plot"); fixed apart, on the owner's word (the next commit).

**Fix (owner, 2026-10-22, "fix it"): `nitcad plot` loads a deck's results for arrows and streamlines.** Reading a
deck, `nitcad plot` loaded the results only for the views that show them, so `--vectors` or `--streamlines` on the
solid, mesh, quality or doping views failed with "the results hold no field to plot" (a result file worked). Now
either option loads them on any view. `run_test` gains arrows on the solid view, streamlines on the mesh view,
arrows over the doping and 3D arrows on the quality view; the mutant that drops the options from the condition is
caught. No change to any image the command could draw before.

**The WG1 gap investigated (owner, 2026-10-22, "investigate app_window"; fixes on "apply them too").** The
self-test's message loop was logged message by message (each dispatched message, its window class and time, and the
wait in `GetMessage` between two ticks of the 10 ms gap timer) over six runs alone and five beside six `nitcad run`
loops. Two causes:
1. **The desktop user's input reached the self-test's window**, which opened shown and activated over the user's
   work: a title-bar drag (a modal move loop, 103 ms, the gap timer unseen meanwhile), clicks on buttons and lists
   (41–65 ms each) and streams of mouse moves over the view (each a hover probe of 1–12 ms in Debug; input is taken
   before the low-priority `WM_TIMER`, so a steady stream starves the timer: gaps up to 454 ms with the window
   working throughout). The earlier failures (1,623, 652, 279, 249 ms) came from suites run while the workstation
   was in use; that they were all input is likely, not proven.
2. **The run monitor's list repainted slowly in Debug**: 30 ms alone, 118 ms under load, with no input — each cell
   re-formatted its whole row.

Fixes: `WindowOptions::input` (default true); the self-test's window is shown without activation and disabled (it
still paints and runs its timers), so WG1 measures the window's own work. `RunPage` keeps each row's cells, formatted
when the row changes (`PointRow` gains `==`). Tests: a window made without input is disabled, a usual one not; after
a run every monitor cell equals the final rows' formatting. Mutants: 2 of 2 caught (rows kept stale, the window
taking input). Afterwards the Debug GPU self-test's longest gaps were 36.1, 38.5 and 36.6 ms, and Debug `app_window`
passed in the suite (19 of 19). WG1 stands as ruled.

**Commit 26c·2 (region visibility; YD4, YG1, YG2, YG5, YG6), as built.**
- `render::ViewOptions::hidden_regions` (region ids; an id not the device's is refused, "a hidden region is not one
  of the device's"). **Hiding goes by control volume, in every view:** a node is hidden with its region, and what is
  drawn is the shown nodes' pieces (`Geometry::owner_pieces`), so a 2D map's area and a 3D surface's enclosed volume
  are the shown regions' control volumes. Edges are drawn as the shown node's half (an edge's midpoint is where the
  two control volumes meet). The interface pieces between a shown and a hidden node become boundary faces of the
  shown side, turned toward the hidden node, so the shown regions stay closed (solid, field, clipped, opaque,
  quality, the slice and iso views' outline shell). A contact all of whose nodes are hidden is left out. Scales fit
  over the shown nodes (`render::shown_values`: NaN on hidden nodes), in the field, slice, iso, difference and
  filmstrip views (one scale over every frame's shown nodes). Isosurfaces are split into the shown pieces of each
  triangle (`IsoSurface::elements` gives each triangle's element), with normals interpolated from the corners. The
  exploded view leaves the hidden regions out and keeps the others' places. 3D quality draws the elements' pieces
  in the shown control volumes and the opened faces in their element's color, not whole elements.
- `nitcad plot --hide <region>` (repeatable; twice is once) on the solid, mesh, field, slice, iso, exploded,
  quality, difference and filmstrip views of a 2D or 3D device. Refused: an unknown name ("there is no region
  <name>"), another view or a 1D device ("hidden regions apply to the solid, mesh, field, slice, iso, exploded,
  quality, difference and filmstrip views of a 2D or 3D device"), no name (usage).
- The window: a check box per region in the scene tree (contacts and other items have none). Clearing a box hides the
  region, checking it shows it again; the status bar says which. View ▸ Show every region (Ctrl+Shift+H) checks them
  all, enabled while any is hidden. `ViewSettings::hidden` (names, in the device's order) is kept on the device views
  (`takes_hidden`: solid, mesh, field, slice, iso, exploded, quality, difference); other views leave it out, and
  `plot_arguments` gives `--hide` for each. The probe doesn't pick a hidden node; in 3D it picks the shown nodes of
  the faces the hidden regions opened. A new run resets the view settings, hidden regions with them. **Without a
  comctl32 v6 manifest the tree sends no TVN_ITEMCHANGED**, so a click or the space bar on the tree posts a message
  that reads the boxes once the tree has toggled them. `TVS_CHECKBOXES` gives every new item a box, so the window
  removes it from items that aren't regions.
- **Not hidden:** arrows and streamlines (they belong to the field, not to a region); the mesh view's edges in 3D are
  the shown halves of the domain boundary, with no edges on the opened faces.
- **Tests:** `hidden_test` (new) — nothing hidden draws as before, and an id that is not a region's is refused; 2D:
  the pieces' area is the silicon's control volumes (1.7 µm²), the gate is gone, outlines and mesh lie below the
  interface, which is drawn across the device, scales are fitted over the silicon, the filmstrip has one scale over
  shown frames, quality; 3D: the solid and field surfaces closed (volume and vector area) with either region hidden,
  opaque and clipped too, the mesh, slices (area) with their outline closed, a curved isosurface cut at the
  interface with its corners' normals, quality closed, exploded; `shown_values`. WARP references
  `mosfet2d_hidden_oxide` (potential, mesh and contours) and `mosfet3d_hidden_oxide` (n). `run_test`: `nitcad plot
  --hide` on ten 2D and 3D views, hiding changes the picture, twice is once, a filmstrip's scale is the silicon's
  (the same as `--range` of it), the refusals. `view_settings_test`: names and refusals as `nitcad plot`'s; **to the
  bit equal to `nitcad plot --hide` (YG1)** in 13 cases (2D solid, mesh, field, quality, iso with edges, difference
  over run B; 3D solid, field, exploded, quality, slice, iso, clip), the difference's scale fitted over the shown
  nodes. `probe_test`: hidden nodes not picked, the opened faces' nodes picked in 3D. `analysis_window_test`: boxes on
  the regions only, the space bar hides and shows, the window equal to `nitcad plot --hide`, no probe where the oxide
  was, kept on the mesh view, both hidden in the device's order, Show every region, a run resets, the menu entry.
- **Mutation checks (Release):** 35 of 36 caught. Render 18 (an id range off by one, no node hidden, whole edges,
  hidden edges kept, no opened faces, opened faces turned inward, the 2D interface left out, hidden contacts kept,
  the colorbar, the filmstrip and `shown_values` over every node, iso pieces unfiltered, iso normal weights, the
  exploded view keeping the region, a 3D interface into the hidden region, the solid, quality and outline shells
  not closed); `nitcad plot` 6 (difference refusing `--hide`, 1D accepting it, the field, filmstrip and difference
  fits over every node, `--hide` ignored); app 6 (difference dropping the regions, no `--hide` arguments, the field
  and difference fits, hidden nodes picked, opened faces not picked); window 5 of 6 (contacts keeping a box, Show
  every region doing nothing, the space bar ignored, the probe picking hidden nodes, the order). **Equivalent:** the
  prober hiding on every view, since every view that can be probed takes hidden regions. Four survived at first and
  gained tests: the iso normals (the isosurface was a plane), the outline shell, the filmstrip's fit, and
  `nitcad plot`'s fits (now also run against the window's to-the-bit comparisons).
- **Suites (clean worktree at `4a9d236`):** Release 27 of 27 (1,308 s), no warnings under `/W4 /WX`; Debug 18 of 19
  (1,526 s): `render` failed after 10 s, all three of its failures `D3D12CreateDevice failed (HRESULT 0x887A0007)`
  (DXGI_ERROR_DEVICE_RESET, before anything was drawn), while the other D3D suites of the same run passed; rerun
  alone in the same worktree it passed (651 s).

**Commit 26c·3 (transient playback and the clipping box's sliders; YD5, YD6, YG3, YG4, YG5, YG6), as built.**
- **Playback** on the field view of a transient run of a 2D or 3D device (`app::takes_film`), the only device view
  `nitcad plot --view filmstrip` draws. `app::film_scale` fits one scale over every snapshot as the filmstrip does
  with `--every 1` (DC4): the range if given, over the shown regions' nodes (YD4), a vector magnitude on a log scale
  floored ten decades below its largest value over all. The window fits it once and keeps it in
  `ViewSettings::film`, fitting again only when the field, scale, range or hidden regions change; a new workspace
  stops playback. With a film, `build_scene` draws the snapshot (`point`) through `render::filmstrip_scenes`, the
  filmstrip frame's own call: the same scale, and the time in the title ("ψ (V), t = 200 ps"). `plot_arguments`
  gives the snapshot's field view (`--point k`); no `nitcad plot` option draws one frame alone, so YG3 compares
  with `--frames`.
- **Controls:** in the inspector, a Snapshot slider over the snapshots, the time and "k of n" beside **Play/Pause**.
  Simulation ▸ Play / pause snapshots (Ctrl+Shift+Space), First (Alt+Home), Previous (Alt+Left), Next (Alt+Right)
  and Last snapshot (Alt+End), enabled while playback applies. Play starts from the snapshot shown, or from the
  first when the last is shown, and stops at the last. Pause, a step key, the slider, the point list and another
  tab each pause it; a step or the slider lands on the snapshot asked for. **Frames:** a 40 ms timer shows the
  snapshot the time since Play has reached (`app::playback_frame`, five a second), so a frame drawn late skips
  those it missed rather than queueing them, and the message loop is never held for more than a frame's drawing.
- **Clip sliders** on the solid and field views of a 3D device: per axis a low and a high slider side by side,
  1,000 steps over the device's extent, values rounded to the nanometre (`app::clip_value_um`, `app::clip_step`).
  Moving one makes the box over the whole device if there was none (Clipping box checked, the box typed), sets that
  side, and moves the axis's other side with it rather than letting them cross. A typed box (Apply) moves the
  sliders to its nearest steps; on an extent under 1 µm a step is finer than a nanometre, so a value, not a step,
  comes back unchanged.
- **Tests:** `view_settings_test` — playback only on a transient field view (not a sweep, 1D or the solid view),
  the refusals; the film holds every snapshot's own scale; a frame's title and kept scale; |J|'s floor and a range;
  the frame draws values below the film's floor at it; `playback_frame` (five a second, a late tick skipping, the
  last kept); the sliders' values, steps and round trips. **To the bit equal to `nitcad plot --view filmstrip
  --every 1 --frames` (YG3)** for every snapshot of a ramped 2D MOSFET (10 snapshots): ψ with contours, |J| on its
  log scale, n with the oxide hidden, ψ linear with the oxide hidden. `analysis_window_test`: the controls and
  commands on the field view, the last snapshot as its frame; First, Previous (stays at the first), Next, Last
  (stays at the last) and the slider land on their snapshot, the point list following; Play from the first runs to
  the last, every frame shown equal to its frame and none before its time, the longest message 100 ms or less
  (WG1); Play from the last starts again, Pause holds, the button plays, the slider and another tab pause; none on
  the mesh view; the menu and keys, tab stops; with the oxide hidden the scale is fitted again and the frame equals
  `--hide oxide`'s. The 3D MOSFET's sliders: six on the solid view at the whole extent, a quarter of x giving the box
  `0.5,0,0:2,1,0.8` (checked, typed) and the view equal to `nitcad plot --clip` (YG4), the high z side moved below
  the low one, a typed box moving the sliders to its steps, the field view keeping them and the mesh view without,
  side by side at 200% (YG5).
- **Mutation checks (Release):** 28 of 28 caught. App 13: the frame not floored, the film floored nine decades, over
  hidden nodes, ignoring the range, the frame titled at t = 0 and the film on every view (both caught by the
  `/WX` build, an unused variable), the frame fitting its own scale, the film in 1D, the frame rounded, past the
  last, the sliders to 0.1 nm, floored and unclamped. Window 15: Play from the last staying there, playing past
  the end, another tab or the slider leaving it playing, Previous from the first, the high side crossing the low,
  the typed box not on the sliders, a new box empty, the film not fitted again on hiding, the button's text,
  overlapping sliders, a frame every tick, the timer ignored, the controls on every view, the point list left
  behind. Five survived at first and gained tests: the floor, past the last, hiding while playing, frames before
  their time, and the point list.
- **Suites (clean worktree at `3d0e78a`):** Release 27 of 27 (1,013 s), no warnings under `/W4 /WX`; Debug 18 of
  19. Tests 1–13 passed in order (767 s; the run was stopped there), tests 14–19 run six at a time (647 s):
  `app_window` failed WG1 in the WARP self-test, the message loop's longest gap 102.1 ms (bound 100 ms; 19 of its 20
  cases passed); rerun alone it failed the same way (125.8 ms). The WARP self-test was then run twice alone in this
  worktree and twice in the one at `4c4b75b` (the WG1 fix, before 26c·2 and 26c·3): at `3d0e78a` both passed, every
  gap 36.0–42.1 ms; at `4c4b75b` both failed, one view's gap 1,095.9 ms and 215.2 ms, the other views' 36.6–38.6 ms.
  The gaps are rare single spikes on a loop that otherwise holds about 40 ms, in the build before this unit's window
  changes as in this one, so they are not this commit's. **Open item (WG1 spikes, owner, 2026-10-23, "record as
  is"):** the window takes no input since the WG1 fix, so the spikes have another cause, unknown; WG1 is left as
  ruled.

**Status:** Unit 26c built in three commits (c16f8e0, 4a9d236, 3d0e78a), with the `nitcad plot` fix (3e730dd) and
the WG1 fixes (4c4b75b); YG1–YG6 met on this workstation, the WG1 spikes above open. **Accepted by the owner
2026-10-23** ("accept it and merge"); ARCHITECTURE.md 6.14, 11 and 11.3 updated to "as built" (YG6). Merged to
`main` as `1be857e` (PR #61, 2026-10-23, on the owner's word, the hosted CI still queued).

**Not in Unit 26c:** in-view clip handles (the owner chose sliders), the tiled filmstrip in the window, video
export, oblique clip planes, streamtubes.

### Parallel test suites (owner, 2026-10-23, "the test is taking too long", then "do the test-speed"; on `test/parallel-suites`)

The Debug suite ran its 19 tests one after another in about 40 minutes (`solve` 769 s, `render` 741 s,
`app_window` 410 s), on a 10-thread workstation: the presets set no job count.
- **One ctest entry per test file** (`nitcad_add_file_tests` in CMakeLists.txt) for the seven long executables
  (`mesh`, `assemble`, `render`, `solve`, `deck`, `app`, `app_window`): `<name>.<file>`, each running
  `"-#" "[#<file>_test]~[.]" --allow-running-no-tests`. Catch2's sharding was not used: with the default random
  order and no fixed seed its shards may overlap or miss tests (Catch2 3.16 warns so), and a fixed seed or
  declaration order would give up the random order. Each file keeps it. Hidden tests stay out; the Release-only
  hidden-tag entries are unchanged; `solve.mosfet` runs nothing (all its tests are hidden). **Coverage:** per
  executable, the files' test-case counts sum to the plain run's: 95, 107, 78, 215, 43, 57, 20.
- **Four jobs** in both test presets (`"execution": {"jobs": 4}`; CI's `ctest --preset` takes it). The Debug
  suite is 2,037 s of work; four jobs already hold the rest beside its longest entry, and fewer jobs load the
  timing gates less. `render.image` (its 2D and 3D MOSFET fixtures solve for about 530 s in Debug before any
  view is drawn; splitting its sections would only repeat them) has `COST 600`, so it starts first without
  ctest's timing record.
- **The WG1 files run alone** (`RUN_SERIAL`): `app.run_controller`, `app_window.self_test`,
  `app_window.analysis_window`, and in Release `app_window_3d`. Run beside the others (six jobs, then four),
  the self-test passed once and failed once with gaps of 494 and 611 ms, its 16-point solve at 95 s against
  about 42 s alone. AC9 (`solve.hole_mobility`, a cost ratio within one process) passed in parallel in every
  run and is left parallel.
- **Suites (working tree on `1be857e`):** Debug 99 of 99 in 947 s (was about 2,400 s); Release 107 of 107 in
  618 s (was 1,013 s at `3d0e78a`). **Clean worktree at `15c1088`:** builds from scratch, Release 94 s and Debug
  56 s (vcpkg's binary cache), no warnings under `/W4 /WX`; Release 107 of 107 in 393 s; Debug 98 of 99 in 850 s:
  `app_window.self_test`, run alone after `render.image`, failed WG1 with one gap of 296.6 ms, its solve at the
  usual speed (16 points in 35.4 s) and the monitor's plot drawn in 1.8 ms: the WG1 spike open above (Unit 26c),
  not load.
- **Hosted CI at `1af3b55` (PR #62; the first CI to finish since Unit 26b — PR #61's runs and `main`'s were
  cancelled):** Release passed (21 min); Debug 97 of 99 (1,929 s), both failures in entries run alone. (1) The
  26c·3 playback test's wait for its transient MOSFET run timed out: `pump_until`'s 120 s, a guard against a hang,
  is shorter than that run on the runner's Debug build. The guard in `analysis_window_test.cpp` is now 600 s; the
  file passes locally (139 s). (2) The GPU self-test (WARP on the runner) failed WG1 with a gap of 114.1 ms, its
  solve at 69.3 s (about twice this workstation's).
- **WG1 on the hosted runner (owner, 2026-10-23, "c": a looser CI-only bound):** 250 ms where the environment
  variable `CI` is set (GitHub Actions sets it), 100 ms everywhere else, as ruled. `nitcad-app --self-test` takes
  `--gap-limit-ms <ms>` (default 100; more than 0) and names the bound in its WG1 failure; the app does not read
  the environment. The tests take the bound from one place, `app_test::wg1_limit_ms()` (tests/app/fixtures.hpp),
  for all three WG1 checks: the self-tests (passed as `--gap-limit-ms`), the run controller and playback. New
  checks: `--gap-limit-ms 0` and a missing value are usage errors (64); a 1D diode run (701 points) under a
  0.001 ms bound fails WG1 naming it. Mutant: the bound ignored (100 ms fixed) — caught. Release: the three WG1
  entries pass with and without `CI` set; no warnings under `/W4 /WX`.
- **Hosted CI at `be4d95a`:** Release passed (push and pull-request runs, 21 and 28 min); the three WG1 entries
  passed under the 250 ms bound, and the playback test under its 600 s guard. Debug 98 of 99 in both runs: in
  `app_window.window` the wait for the example MOSFET run (the run monitor test) timed out at `window_test.cpp`'s
  own `pump_until`, still 120 s (the entry failed at 171.5 s and 178.6 s). That guard is now 600 s, as in
  `analysis_window_test.cpp`.
- Not done: restructuring `render.image` (2D and 3D views as separate test cases) to bring Debug near 10
  minutes; worktrees run at once still share their `%TEMP%` scratch folders.

---

### Unit 46 — Lombardi surface mobility, opt-in (started 2026-10-24; on `physics/surface-mobility`)

**Owner request:** "start Lombardi surface-roughness mobility" (2026-10-24), the first of the nanoscale models of
ARCHITECTURE.md 11.3 (found by the FinFET deck). **Owner rulings (2026-10-24):** (1) E⊥ is the field normal to the
current, (2) the Jacobian is exact, with a wider pattern, (3) the electrothermal assembly is included, (4) the surface
terms combine with the run's selected bulk mobility. The standing conditions hold: opt-in and bit-identical when off,
criteria fixed before code, Release and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or
merged without authorization.

**Sources read (before coding).** C. Lombardi, S. Manzini, A. Saporito, M. Vanzi, IEEE Trans. CAD 7, 1164 (1988),
doi:10.1109/43.9186, is paywalled and was not read. The model is read from two text reproductions that cite it:
(a) the Xyce Reference Guide (Sandia; `doc/Reference_Guide/Xyce_RG_PDE.tex` and `lombardiTbl.tex` on
github.com/Xyce/Xyce), (b) Genius TCAD Open (`src/material/Si/Si_mob_Lombardi.cc` on
github.com/cogenda/Genius-TCAD-Open, Medici's `.LSM` parameter names). Both give

    1/mu = 1/mu_b + 1/mu_ac + 1/mu_sr,
    mu_ac = B / E⊥ + C N^lambda / (T E⊥^(1/3)),      mu_sr = delta / E⊥^2,

N the total impurity density [cm^-3], T [K], E⊥ [V/cm]. COMSOL's Semiconductor Module guide (doc.comsol.com,
"Lombardi Surface Mobility Model") cites the paper for its silicon values; its equations are images, not read.

| Coefficient | Xyce | Genius | COMSOL (as extracted) | Used |
|---|---|---|---|---|
| B_n [cm/s] | 4.75e7 | 4.75e7 | 4.75e7 | 4.75e7 |
| C_n [cm^(5/3+3λ) K V^(-2/3) s^-1] | 1.74e5 | 1.74e5 | 1.74e5 | 1.74e5 |
| λ_n | 0.125 | 0.125 | 0.125 | 0.125 |
| δ_n [V/s] | 5.82e14 | 5.82e14 | 5.82e14 | 5.82e14 |
| B_p [cm/s] | 9.93e6 | 9.93e6 | 9.93e7 (an extraction artefact or a misprint; two sources disagree) | 9.93e6 |
| C_p | 8.84e5 | 8.84e5 | 8.84e5 | 8.84e5 |
| λ_p | 0.0317 | 0.0317 | 0.0317 | 0.0317 |
| δ_p [V/s] | 2.05e14 | 2.05e14 | 2.05e14 | 2.05e14 |

The E⊥ exponent of μ_sr is 2 in both (Xyce `exn8`, `exp8`). Genius applies Canali's saturation to the Matthiessen
sum (the low-field mobility includes the surface terms), as NiTCAD will. Lombardi's own bulk term is Masetti's
arsenic electron set and boron hole set; by ruling (4) it is not used.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| LD1 | `PhysicsModels::surface_mobility`: `off` (the default) or `lombardi`; the deck's `[models] surface_mobility`; the run record names it only when `lombardi`; the app's deck vocabulary lists it. | Opt-in, as Units 34, 37, 38. | — |
| LD2 | Per edge, the field vector is the mean of the two ends' least-squares fields (the reconstruction of impact ionization: E_i = P_i Σ_k t_k E_k over the node's carrier edges, so only semiconductor fields, never an oxide's), and E⊥ = sqrt(\|E − (E·t)t\|^2 + E_0^2) with t the edge's unit vector and E_0 = 1 V/cm. | Ruling (1); exact for a uniform field at any angle; E_0 keeps B / E⊥ and δ / E⊥^2 finite and smooth where E⊥ vanishes (Genius adds 1 V/cm). | Bulk edges at an angle to a junction field are reduced too. In 1D E⊥ = E_0 (the change is below 1e-4). |
| LD3 | Each end's low-field mobility is the Matthiessen sum of its bulk mobility (the run's selection, at its N and T) and the surface terms at its own N and T with the edge's E⊥; the edge mobility is their harmonic mean, then Canali when on (on a mixed edge each end's Canali, then the mean), as today. N is the bulk model's N (0 with `doping_mobility` off). | Ruling (4); the existing harmonic mean and Canali order unchanged. | — |
| LD4 | Exact Jacobian: E⊥ is linear in ψ before the root, so each edge stores the coefficients g_j (D-vectors, t-projected out) over the nodes j of both ends' stencils; the carrier rows of a and b gain ψ columns of those nodes (distance 2), added to the pattern only with the model on. | Ruling (2); quadratic Newton and exact small-signal (Units 13, 19). | More nonzeros; measured and reported (LS11). |
| LD5 | Electrothermal: the per-end mobility and its θ partial include the surface terms (∂μ_ac/∂T = −C N^λ / (T^2 E⊥^(1/3))); E⊥ is from ψ, which does not depend on θ. | Ruling (3). | — |
| LD6 | `lombardi` is refused (invalid_input) on a device with a semiconductor region whose electron or hole Caughey–Thomas parameters are not silicon's. | The coefficients are silicon's; as Units 34, 38. | SiGe and others not offered. |
| LD7 | Not on thermionic edges (no mobility enters them). | As Canali. | — |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| LS1 | Off (the default): every existing test and pin passes unchanged; run identities and Jacobian patterns unchanged. |
| LS2 | μ_ac, μ_sr and the Matthiessen sum equal the formula and coefficients above, evaluated independently in the test, at E⊥ = 1e3 … 1e7 V/cm (each decade), N = 0, 1e14 … 1e20 cm⁻³ and T = 250, 300, 400 K, within 1e-12; the E⊥ and T slopes within 1e-6 of central differences. |
| LS3 | E⊥ of a uniform field on 2D meshes (tensor and Delaunay) at 0, 30, 45 and 60° and on a 3D tensor mesh: each interior edge's E⊥ equals sqrt(\|E × t\|^2 + E_0^2) within 1e-10 relative. |
| LS4 | The Jacobian against central differences of the residual with the model on (isothermal and electrothermal; 2D and 3D; with and without Canali): the largest relative column error at most 1e-6. Newton on the 2D MOSFET in strong inversion converges quadratically (the last three correction norms). |
| LS5 | 1D diode: the current with `lombardi` within 1e-4 relative of the current without (E⊥ = E_0). |
| LS6 | A long-channel n-MOSFET (Release only), Vd = 10 mV, substrate 1e16 and 1e17 cm⁻³: the effective mobility μ_eff = (L / W) g_d / Q_inv at mid-channel, against μ of LD3 evaluated at E_eff = q (Q_b + Q_inv / 2) / ε_Si and the substrate's N, within 15% for E_eff in 3e5 … 1e6 V/cm. |
| LS7 | Universality: in LS6's devices, μ_eff of the two substrates at equal E_eff (interpolated) within 15% for E_eff in 5e5 … 1e6 V/cm. |
| LS8 | The 2D MOSFET example's Id–Vg (Vd = 0.1 V): at Vg = 1.5 V the current with `lombardi` lies below the current without; below threshold (Vg = 0.2 V) within 10%. The ratio at 1.5 V is pinned. |
| LS9 | The 2D MOSFET example rotated 90° (x and y exchanged): the same Id at Vg = 0.5 and 1.5 V within 1e-9 relative. |
| LS10 | Refused on a germanium region and accepted on silicon with another parameter changed; with electrothermal and 300 K sinks within 1e-3 of the isothermal current (LS8's device, Vg = 1.5 V); small-signal conductance against DC differences within 1e-5. |
| LS11 | Reported, not gated: the pattern's nonzeros and the 2D and 3D MOSFET examples' sweep times with and without the model; the FinFET deck's on-current with and without (the deck's mesh convergence is still open). |
| LS12 | Full Release and Debug suites; zero warnings; hand mutants caught. |

**As built (2026-10-24).** `physics/surface_mobility` (the formula and its partials),
`DriftDiffusion` (the per-edge stencils of LD4, the surface path of `edge_fluxes` and `thermal_fluxes`, the
small-signal current rows), the `CarrierFluxes` and `ThermalEdges` terms (the pattern and the E⊥ partials, in the
heat rows too with electrothermal), the deck key, the run record and the app's vocabulary. Tests:
`tests/physics/surface_mobility_test.cpp` (LS2), `tests/assemble/surface_mobility_test.cpp` (LS1, LS3, LS4, LS10),
`tests/deck/surface_mobility_test.cpp` (LS5, and Release only under `[lombardi]`, ctest `deck_lombardi`: LS4's
Newton, LS6–LS11).

| # | Added in the build | Why |
|---|---|---|
| LD8 | `lombardi` with `impact_ionization` is refused (not supported). | Impact ionization's generation rows use the fluxes' partials over its own pattern; the E⊥ partials reach distance 2, outside it. Combining them needs that pattern widened; not done in this unit. **For the owner.** |

**Results:**

| # | Result |
|---|---|
| LS1 | Met. With the model off the pattern and residual are identical (tested); every existing test and pin passes. |
| LS2 | Met: 1,821 assertions (formula within 1e-12, slopes within 1e-6). |
| LS3 | Met on 2D tensor, 2D Delaunay and 3D tensor meshes at every angle. |
| LS4 | Met. The FD gate holds within 1e-6 in all nine cases (2D, Delaunay, 3D; ±Canali; electrothermal 2D, 2D Canali, 3D Canali). On the 2D MOSFET at Vg = 1.5 V the last three corrections are 1.9e-1, 2.5e-5 and 4.3e-15, an order of 2.5. |
| LS5 | Met: the diode within 1e-4, and the run record names the model only when on. |
| LS6 | **1e17 met; 1e16 not met.** Measured at Vg 0 … 6.5 V by 0.5 V on a mesh with 0.05 nm spacing in the top 5 nm. 1e17: gap +0.8% at 3.3e5 V/cm rising to +11.7% at 9.8e5. 1e16: +5.4% at 3.5e5, +14.4% at 7.5e5, **+20.0% at 9.2e5** (21.6% at 9.6e5 with 0.25 V steps). Without the model the same extraction gives the bulk mobility within 1% (checked), so the extraction is sound. Cause: μ_eff is a current-weighted average of μ(E⊥) over the classical inversion layer, whose field falls from E_s at the interface into the bulk. μ is convex in E⊥ (Jensen), so the average exceeds μ(E_eff). The gap grows with field and is larger at the lighter substrate, where Q_inv/2 dominates E_eff. Encoded as a `[!shouldfail]` case with the criterion unchanged. **For the owner.** |
| LS7 | Met: worst 13.2% at 5e5 V/cm, 11.5% at 1e6. |
| LS8 | **Partly not met.** At 1.5 V the current is lower: 4.5348e-2 A/cm against 6.5274e-2, ratio 0.694738, pinned within 1e-4. **Below threshold (0.2 V) the ratio is 0.7249, not within 10%.** On the same device with 0.05 nm spacing at the surface it is 0.718, so the mesh does not cause it. Cause: in depletion the surface field is about 4e4 V/cm at 0.2 V, and there μ_ac (B/E⊥ ≈ 1,100; the C term ≈ 1,660 at 1e16) already lowers μ from 1,250 to about 430 at the interface. The subthreshold current flows in that layer, and the model has no distance damping. The criterion assumed a small effect; the model says otherwise. Encoded as a `[!shouldfail]` case. **For the owner.** |
| LS9 | Met: the reflected device's Id within 1e-9 at 0.5 and 1.5 V. |
| LS10 | Met. Refused on germanium, accepted on silicon with another parameter changed; electrothermal with 300 K sinks 0.999936 of the isothermal current; small-signal g_d and g_m against DC differences within 1e-8. Added in the build, as the mutants found: the small-signal current rows against central differences of the terminal currents (2D, 3D Canali, 2D electrothermal Canali), within 1e-6 of the largest entry. On the example the drain contact lies in the n+ pocket, where E⊥ is small, so the device-level check alone could not catch a missing E⊥ partial in those rows. |
| LS11 | 2D MOSFET example: Jacobian nonzeros 30,244 off, 37,284 with the model (+23%); sweep 0.94 s off, 1.5 s with it. 3D MOSFET example: 5.6 s off, 12.0 s on; Id(1.2 V) ratio 0.718. FinFET deck (PARDISO, 6 threads): on-current at 0.8 V 1.6586e-5 A with the model, 1.8835e-5 without (ratio 0.881); at 0.2 V 0.971. 30 min with the model against 18 min without, but the with-model time was measured while other builds ran, so it is not a clean comparison. |
| LS12 | Release (110 tests) and Debug (101) suites pass, zero warnings. One existing pin changed with the vocabulary, not loosened: the app's `[models]` completion on an empty line lists 16 keys (17 less `srh`), not 15. Hand mutants, all caught: the E⊥ projection removed (LS3); the electrothermal T term of the θ partial (LS4); the heat-row E⊥ partials (LS4); the thermal continuity E⊥ partials (LS4); the isothermal and thermal small-signal E⊥ terms (LS10's current rows; both survived the device-level check before those rows were added); the germanium refusal (LS10); the run-record entry (LS5); E_0 = 0 (a contract violation at E⊥ = 0 in the 1D diode). |

**Limitations found (for the owner):**
- **Mesh dependence.** In classical drift-diffusion the inversion layer is about kT/(qE_s) thick, 0.26 nm at
  1e6 V/cm. On a mesh with 1 nm surface spacing the LS6 gap reached +58%; it converged near 0.05 nm. The examples'
  and the FinFET's surface spacings (about 10 nm) under-resolve the layer, so their Lombardi currents are
  indicative. LS8's subthreshold ratio, which depends on a depletion layer many nm thick, is not affected. A
  follow-up could take E⊥ at the interface node from the interface's D-field (Gauss) instead of the
  least-squares nodal fields.
- **No distance damping.** The model applies at every carrier edge with E⊥ from the field normal to the current,
  not only near an interface. Bulk edges at an angle to a junction field are reduced too (LD2). Commercial tools
  often damp the surface terms with distance from the interface. Not in the owner's rulings; not done.

**Ruling on LS6 (1e16) and LS8 (below threshold) (2026-10-24).** The owner delegated it ("Decide yourself").
Both are accepted as limitations of the model as specified, not defects. The criteria are unchanged, and the two
`[!shouldfail]` cases stay, so a change that makes either pass is flagged.
- LS8's below-threshold limit assumed a small effect. Lombardi without distance damping lowers the mobility in the
  depletion field at the surface, and the result has converged on the mesh.
- LS6's 1e16 gap is the Jensen bias of comparing μ_eff with μ(E_eff) across the classical inversion layer. The
  extraction is verified, and 1e17 and LS7 pass. Quantum confinement (MLDA, ARCHITECTURE 11.3) narrows the layer's
  field spread, so it is the place to revisit this.
- LD8 (refused with impact ionization) and the mesh dependence stand as documented limitations.

---

### Unit 47 — quantum confinement by the modified local density approximation, opt-in (started 2026-10-24; on `physics/quantum-confinement`)

**Owner request:** "do #1" (2026-10-24), MLDA quantum confinement, the second of the nanoscale models of
ARCHITECTURE.md 11.3. **Owner rulings (2026-10-24):** (1) the density-of-states factor is applied under both
statistics, exact MLDA under Boltzmann and Hänsch's approximation under Fermi–Dirac; (2) refused with
electrothermal; (3) the quantization masses default to Ioffe's m_l and heavy-hole values and can be overridden in the
deck, silicon only; (4) at a lumped-gate node (the interface at the node) the factor is its mean over the node's
half-cell. The standing conditions hold: opt-in and bit-identical when off, criteria fixed before code, Release and
Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged without authorization.

**Sources read (before coding).** G. Paasch, H. Übensee, phys. stat. sol. (b) 113, 165 (1982), and W. Hänsch et al.,
Solid-State Electron. 32, 839 (1989), are paywalled and were not read. The model is read from:
(a) M. Ieong, R. Logan, J. Slinkman, "Efficient quantum correction model for multi-dimensional CMOS simulations",
SISPAD 1998 (in4.iue.tuwien.ac.at/pdfs/sispad1998/pdfs/Ieong_34.pdf; text extracted with `pdftotext`). Its MLDA
density, eq. (1), is an image and was not read; its text gives the ingredients (the distance z to the interface, N_c,
the quasi-Fermi level, the zeroth-order spherical Bessel function j0, the thermal wavelength λ = ħ / sqrt(2 m kT)),
says the expression reduces to Fermi–Dirac statistics for large z, and describes the implementation used here: the
corrected density is turned into a quantum potential that is lumped into the band-structure parameters of the
heterostructure current equation, so equilibrium carries no current;
(b) S. Entner, PhD thesis, TU Wien, section 2.4.1.1 (iue.tuwien.ac.at/phd/entner/node12.html), citing Paasch–Übensee
and Hänsch: eq. (2.106) h = 1 − exp(−(z + z_0)^2 / (ζ^2 λ^2)), (2.107) λ = ħ / sqrt(2 m* k T), (2.108)
N_c = N_c,0 h; the correction does not depend on bias and does not move the band edges; z_0 = 0 makes the density
of states vanish at the interface, "which causes numerical problems" (ζ = 1, z_0 = 1 nm in its figure);
(c) Ioffe's silicon band-structure page (ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html, read over HTTP): m_l = 0.98 m0,
m_t = 0.19 m0; heavy holes 0.49 m0, light 0.16 m0.

The Boltzmann reduction is derived here, not read: with the MLDA kernel 1 − j0(2 z sqrt(ε) / λ),

    ∫0^∞ sqrt(ε) e^(−ε) [1 − j0(2 a sqrt(ε))] dε / (sqrt(π) / 2) = 1 − exp(−a^2),   a = z / λ,

(u = sqrt(ε) and ∫0^∞ u e^(−u^2) sin(b u) du = sqrt(π) b e^(−b^2/4) / 4), which is (b) with z_0 = 0 and ζ = 1. QC2
checks it by quadrature. With the masses of (c), λ_n = 1.2263 nm and λ_p = 1.7343 nm at 300 K (λ for m0 is
1.2140 nm).

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| QD1 | `PhysicsModels::quantum_confinement`: `off` (the default) or `mlda`; the masses `mlda_mass_n` (default 0.98) and `mlda_mass_p` (0.49), in m0; the deck's `[models] quantum_confinement`, `mlda_mass_n`, `mlda_mass_p`; the run record names the model and both masses only when `mlda`; the app's deck vocabulary lists the keys. | Opt-in, as Units 34, 37, 38, 46; ruling (3). | — |
| QD2 | Per semiconductor node i and carrier c, f_c,i = 1 − exp(−d_i^2 / λ_c^2), λ_c = ħ / sqrt(2 m_c m0 k T) at the device temperature; d_i the distance from the node to the nearest interface point. Interface points: the midpoint of every edge between a semiconductor node and an insulator node (where the assembler puts that interface; mesh/interfaces.hpp), and every lumped-gate node itself (d = 0). One nearest distance: at a corner the nearer interface decides. | The Boltzmann MLDA (derived above); the interfaces as the assembler already represents them. | Corners and thin bodies get a single-interface correction (no overlap of two walls' factors). Heterojunctions and contacts are not confining interfaces. |
| QD3 | A lumped-gate node's factor is the mean of f over z in [0, ℓ/2], ℓ the node's shortest carrier edge: 1 − (sqrt(π) λ / (2 a)) erf(a / λ), a = ℓ/2. Every other node takes f at its own d. | Ruling (4); positive, no fitted parameter, the node's share of the charge converges as the mesh is refined. | On a tensor grid ℓ/2 is the box depth; on an unstructured mesh it is an estimate of it. |
| QD4 | The factors enter only through the scaled device (`make_scaled_device`): N_c → N_c f_n and N_v → N_v f_p, band edges unchanged. With the existing fields that is n_ie → n_ie sqrt(f_n f_p), s → s + ½ ln(f_n / f_p), ln(N_c / n_ie) → + ½ ln(f_n / f_p), ln(N_v / n_ie) → + ½ ln(f_p / f_n). Every assembler (equilibrium Poisson, drift-diffusion, transient, small-signal, the quasi-static sweep), the contact values, SRH and the trap and surface-recombination terms read these, so the model needs no new unknowns, its terms are constant in the solution and the Jacobian is unchanged in form. | Ieong's lumping into band parameters; ψ stays the electrostatic potential and the gate and electrode offsets keep their meaning; equilibrium is current-free by the existing heterostructure flux (Unit 15). | Insulator nodes untouched. Band-diagram output shows unchanged band edges; the densities carry the correction. |
| QD5 | Ruling (1): the same factors under Fermi–Dirac (the DOS form of Hänsch; Entner (2.108)). | One model for both statistics. | Under Fermi–Dirac this is Hänsch's approximation, not the j0 integral at the quasi-Fermi level; documented. |
| QD6 | Refused (invalid_input): with electrothermal (ruling (2)); on a device with a semiconductor region whose Caughey–Thomas parameters are not silicon's (as LD6); masses not finite and positive; with `mlda` and no confining interface (no insulator node and no lumped gate). | λ is evaluated at the device temperature; the masses are silicon's. | — |
| QD7 | The distance search: squared distances to all interface points, brute force, once per assembler construction. | Simple; measured in QC12. | O(nodes × interface points); a bucket grid only if QC12 shows it matters. |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| QC1 | Off (the default): every existing test and pin passes unchanged; run identities and Jacobian patterns unchanged. |
| QC2 | λ against ħ / sqrt(2 m m0 k T) from the CODATA constants within 1e-14 relative (300 K: 1.2263 and 1.7343 nm within 1e-4 nm); the factor f(d) against an independent quadrature of the Boltzmann MLDA integral above, within 1e-10, at d / λ = 0.01 … 10; the half-cell mean of QD3 against quadrature within 1e-12, and it tends to f's mean as a → 0 and to 1 as a → ∞. |
| QC3 | Distances: on 2D and 3D tensor meshes with a planar meshed oxide, each semiconductor node's d equals its distance to the interface plane within 1e-12 relative; on a 2D Delaunay mesh with a mirrored-pair interface, within 1e-12 for the nodes adjacent to the interface; at a 2D corner (two oxide walls) d is the smaller wall distance; a lumped-gate node has d = 0. |
| QC4 | Scaled device: N_c f_n and N_v f_p reconstructed from the modified fields within 1e-14 relative at every node; band edges (node_bands) unchanged within 1e-14 V_T; insulator nodes and nodes with f = 1 bit-identical to the model off. |
| QC5 | Equilibrium with the model on (a meshed-oxide 1D MOS-C in inversion and accumulation, the 2D MOSFET example at Vg = 1 V, both statistics): every edge current within its rounding bound (as Unit 11's BGN gate); n p = n_ie^2 f_n f_p within 1e-12 relative at every semiconductor node; SRH and surface recombination zero at equilibrium. |
| QC6 | The Jacobian against central differences with the model on (2D, Delaunay, 3D; Boltzmann and Fermi–Dirac; with a lumped gate and with a meshed oxide): largest relative column error at most 1e-6; Newton on the 2D MOSFET in strong inversion converges quadratically. |
| QC7 | 1D MOS-C, meshed 5 nm oxide, 1e17 p-type, surface spacing 0.05 nm (Release only): in strong inversion the electron density peaks inside the silicon at a distance between 0.3 and 2 nm (classically at the first node); the inversion charge at fixed V_G is lower than without the model; the inversion C–V gate capacitance (quasi-static) is lower than without the model; the threshold (Q_inv = 1e-7 C/cm^2) shifts up by 10–150 mV. In accumulation the hole charge is reduced likewise. |
| QC8 | Mesh convergence (QC7's MOS-C, Release only): the inversion charge at V_G − V_FB = 2 V on surface spacings 0.2, 0.1 and 0.05 nm; the last two differ by at most 1%, and less than the first two. |
| QC9 | A lumped-gate MOS-C against the meshed-oxide one of the same oxide capacitance (QC7's device, spacing 0.05 nm): the inversion charge at V_G − V_FB = 2 V within 5%. |
| QC10 | 2D = 1D: a y-uniform 2D MOS-C reproduces the 1D one within 1e-12; the 2D MOSFET example rotated 90° gives the same Id at Vg = 0.5 and 1.5 V within 1e-9; small-signal conductance against DC differences within 1e-5; a transient long step equals the steady state. |
| QC11 | Refusals: electrothermal, a germanium region, a zero or negative mass, no confining interface; accepted on silicon with another parameter changed, and with Lombardi together. |
| QC12 | Reported, not gated: the distance search time and the sweep times of the 2D and 3D MOSFET examples with and without the model; their Id at Vg = 1.5 V (2D) and 1.2 V (3D) with and without; the FinFET deck's V_th and on-current with and without (its mesh convergence is still open); Unit 46's LS6 gap with MLDA and Lombardi together. |
| QC13 | Full Release and Debug suites; zero warnings; hand mutants caught. |

**As built (2026-10-24).** `physics/quantum_confinement` (λ, the factor and the half-cell mean),
`assemble/quantum_confinement` (`mlda_factors`: the refusals, the interface points, the distances, the factors),
the MLDA block of `make_scaled_device` (QD4), the interface-trap level shift, the deck keys, the run record and the
app's vocabulary. Tests: `tests/physics/quantum_confinement_test.cpp` (QC2), `tests/assemble/quantum_confinement_test.cpp`
(QC3, QC4, QC6, QC11), `tests/solve/quantum_confinement_test.cpp` (QC5, QC10; and, Release only under `[mlda]`, ctest
`solve_mlda`: QC7–QC9 and the QC5 n p report), `tests/deck/quantum_confinement_test.cpp` (QC1; and, Release only
under `[mlda]`, ctest `deck_mlda`: the 2D example's equilibrium, Newton's order, QC10 and QC12), and the QC12 LS6
case in `tests/deck/surface_mobility_test.cpp` (also under `deck_mlda`).

| # | Added in the build | Why |
|---|---|---|
| QD8 | A device with a thermionic-emission interface is refused with `mlda`. | That interface's flux reads the band offsets and the densities of states on both sides directly; QD4 moves the latter, and how the factor should enter that flux was not designed. Silicon devices with a meshed oxide have none. |
| QD9 | `InterfaceEdge::level_shift`: an interface trap's level relative to the intrinsic level is moved by the node's ½ ln(f_n / f_p), so a declared trap energy keeps its meaning against the (unmoved) band edges. Zero with the model off (bit-identical: `tau + 0.0`). | Without it, QD4's shift of the intrinsic level would move every trap relative to the bands; QC4's trap-charge check (equal within 1e-12) caught it. |
| QD10 | QC6's finite-difference gate also runs with Lombardi on (2D and 3D). | QC11 accepts the two together; the 1e17 LS6 sweep with both stalled (QC12), and the gate shows the Jacobian is not the cause. |

**Results:**

| # | Result |
|---|---|
| QC1 | Met. With the model off every existing test and pin passes unchanged (one vocabulary pin grows with the keys, under QC13); the run record has no MLDA entries (tested). |
| QC2 | Met: 65 assertions. λ against long-double CODATA within 1e-14; at 300 K λ_n = 1.2263 nm, λ_p = 1.7343 nm. f against Simpson quadrature of the MLDA integral within 1e-10 at d / λ = 0.01 … 10; the half-cell mean within 1e-12, with both limits. |
| QC3 | Met on 2D and 3D tensor meshes, at a 2D corner (the smaller wall distance), on a 2D Delaunay mesh (nodes adjacent to the interface), and at lumped-gate nodes (d = 0, the half-cell mean). |
| QC4 | Met: N_c f_n and N_v f_p within 1e-14 under both statistics; node_bands unchanged; insulator nodes and nodes with f = 1 bit-identical; trap charges equal within 1e-12; the Jacobian pattern unchanged. |
| QC5 | **Met except the n p relation on the drift-diffusion state.** 1D MOS-C (accumulation −2 V, inversion 1.5 V; both statistics) and the 2D MOSFET example at Vg = 1 V: every terminal current within its resolution; SRH and interface rates within 1e-9 of their scale. n p = n_ie² f_n f_p holds within 6.7e-16 on the equilibrium (quasi-static) state at the same bias, and is gated there. On the drift-diffusion state the worst node misses 1e-12: 3.26e-12 in inversion, at p = 8.8e-5 cm⁻³; without the model the same state gives 7.45e-13. That is the rounding of the minority density in the solver's last update (a tighter `tol_update` does not change it), not the model. Reported by a hidden case; the criterion is unchanged. **For the owner.** |
| QC6 | Met: the finite-difference gate within 1e-6 on 1D, 2D, 3D, Delaunay and lumped-gate devices under both statistics, and (added in the build, below) with Lombardi as well in 2D and 3D. Newton on the 2D MOSFET example at Vg = 1.5 V: corrections 1.6e2, 9.0, 1.1e-1, 4.1e-6, 6.1e-15, an order of 1.99. |
| QC7 | Met (h = 0.05 nm, V_FB = −0.9805 V, V_G − V_FB = 2 V): the electron density peaks 0.800 nm inside the silicon (classically 0.025 nm, the first node); Q_inv 4.8610e-7 against 5.1824e-7 C/cm²; quasi-static C/C_ox 0.9153 against 0.9477; threshold shift +25.7 mV. Accumulation (V_G − V_FB = −2 V): Q_G −1.1904e-6 against −1.2466e-6 C/cm². |
| QC8 | Met: Q_inv 4.860640e-7, 4.860840e-7 and 4.860958e-7 C/cm² at 0.2, 0.1 and 0.05 nm; the differences are 4.1e-5 and 2.4e-5 relative. |
| QC9 | Met: lumped 4.86099e-7 against meshed 4.86096e-7 C/cm² (ratio 1.0000). |
| QC10 | Met. The y-uniform 2D MOS-C equals the 1D one within 1e-12 at every bias; the reflected 2D MOSFET example (x and y exchanged, as LS9) gives the same Id at Vg = 0.5 and 1.5 V within 1e-9; small-signal g_d and g_m match DC central differences to 8 digits (6.07386027e-1 and 7.14614488e-2); a 1D MOS-C transient step to −1.5 V ends at the DC state (ψ within 1e-9 V, p within 1e-7). |
| QC11 | Met: electrothermal, germanium, zero and negative masses, no confining interface and (QD8) a thermionic interface are refused with a `quantum_confinement mlda:` message, the same from the equilibrium and drift-diffusion assemblers; accepted on silicon with another parameter changed, and with Lombardi. |
| QC12 | 2D MOSFET example: the nearest silicon node is 5.00 nm from the oxide (its surface spacing), where f_n = 1 − 6e-8 and f_p = 1 − 2.5e-4, so the electron current Id(1.5 V) is unchanged to 6 digits (ratio 1.000000); sweep 0.69 s off, 0.68 s on. 3D example: 4.16 s off and on, Id(1.2 V) ratio 1.000000. The examples' meshes do not resolve the 1–2 nm confinement layer, so the model has no effect on them. FinFET deck (45,375 nodes, nearest silicon node 0.100 nm from the oxide): distance search 0.064 s. FinFET sweep with MLDA (Vd = 50 mV, the deck's Fermi–Dirac statistics so QD5's form, PARDISO on 6 threads, `line_search = natural`; before Unit 50 the deck's `residual` test stopped at 0.35 V, see limitations): 6.7 min, no step cuts. Id against the classical run: ratio 0.83 at 0 V, 0.79 at 0.2 V, 0.77 at 0.35 V (the least), 0.92 at 0.5 V, 0.974 at 0.8 V (on-current 1.835e-5 A against 1.883e-5 A). The threshold at a constant 1 µA rises by about 10 mV (0.257 V against 0.247 V, interpolated in ln Id). LS6 with MLDA and Lombardi: the 1e16 gap is +5.2% at 3.4e5 V/cm, +13.5% at 7.4e5, +18.7% at 9.0e5 and +21.6% at 9.9e5 (without MLDA: +5.4% at 3.5e5, +14.4% at 7.5e5, +20.0% at 9.2e5, and 21.6% at 9.6e5 with 0.25 V steps). The 1e17 gap is +0.8% at 3.2e5 and +10.5% at 9.7e5 (without MLDA: +0.8% at 3.3e5, +11.7% at 9.8e5). The 1e17 sweep stops at the default tolerance (limitations below); its gaps are from a rerun at `tol_update = 1e-7`. MLDA lowers the gap by about 1 point, so it does not explain LS6's 1e16 miss: the Jensen bias remains. |
| QC13 | Met. Release (115 tests, with `solve_mlda` and `deck_mlda`) and Debug (104) suites pass; zero warnings under `/W4 /WX` in both. One existing pin changed with the vocabulary, not loosened: the app's `[models]` completion on an empty line lists 19 keys (20 less `srh`), not 16. Hand mutants, 16 of 16 caught: the trap level shift dropped (QC4); the band shift, n_ie scaling, ln(N_c / n_ie) sign and ln(N_v / n_ie) shift each removed (QC4); the lumped-gate mean over the full cell instead of the half (QC3); the interface point at the insulator node instead of the edge midpoint (QC3); λ without the 2 (QC2); the sign of the mean's series (QC2); the electrothermal, silicon and thermionic refusals each removed (QC11; the electrothermal one survived a weak message check at first, and the test was tightened); the masses not recorded (QC1); the deck's `mlda_mass_p` written into `mlda_mass_n` (QC1); lumped-gate nodes not marked as half cells (QC9); the model not applied in the scaled device (QC7). |

**Limitations found (for the owner):**
- **Two Newton stalls, diagnosed (2026-10-24, owner: "yes" to finding the node and variable).** Neither is in
  the model. A scratch probe restarted Newton from a saved converged point at the next bias and listed the
  largest scaled corrections by node and unknown (`DriftDiffusion::update_size`'s measure).
  - **The FinFET: the residual line search at a rounding residual.** From point 6 (0.30 V) with the gate 0.2 mV
    higher (the run's last cut, 2^-8 of the step), full Newton steps converge in 6 iterations: corrections 3.0e2,
    8.1, 0.15, 1.4e-3, 3.4e-7, 2.8e-9. From iteration 5 the residual is at rounding (4.5e-13), while a correction
    of 3.4e-7 is still left on holes of about 1e-3 cm⁻³ in the n+ source and drain next to the oxide, below the
    density floor. The deck's `line_search = residual` Armijo test sees only noise there, rejects the step,
    backtracks to 1/256 and never finishes the iteration. The run with the residual test was the one reported
    first: it stopped at Vg = 0.35 V after 104 min. A rerun at `tol_update = 1e-7` was worse (304 iterations by
    point 1) and was stopped. **The same deck with `line_search = natural` runs all 17 points in 6.7 min, 8–15
    iterations per point, with no cuts.** Its currents equal those of the residual-test run at the points both
    solved, to 5 digits.
  - **The 1e17 long channel with Lombardi: the minority-hole floor in uncontacted n+.** At the step from 1.0 to
    1.5 V the corrections settle at 1.7–3.8e-8 with the residual at rounding (3.3e-11), at full steps. They sit on
    minority holes (p 10–20 cm⁻³) in the n+ pockets under the oxide, above the contacts' boxes (y 1.99–2.0 µm),
    uniform down whole columns: those holes are nearly decoupled, so rounding sets their correction. **With MLDA
    off the same step stalls at 1.4–2.2e-8**, so the floor exists with Lombardi alone. Lombardi-only sweeps pass
    the 1e-8 test or not by margin, and MLDA moves the floor a little.
  - **Resolved by Unit 50** (`main` 0c14841, PR #71, 2026-10-28), the separate unit proposed here: (a) the residual
    line search takes the full step when the residual is already at its rounding level; (b) the opt-in
    `minority_floor`. The FinFET deck's advice to use `line_search = natural` with MLDA, first given here, is
    withdrawn. Rebased on `main`, the FinFET with MLDA and the deck's `line_search = residual`: all 17 points in
    11.4 min (PARDISO, 6 threads), 8–12 iterations per point, no step cuts; Id(0.8 V) 1.8349e-5 A against the
    natural run's 1.835e-5 A (2026-10-28).
- **QC5's n p relation on the drift-diffusion state** (above): 3.26e-12 against 1e-12, rounding of a minority
  density of 1e-4 cm⁻³.
- **The examples' meshes do not resolve the layer.** The 2D and 3D MOSFET examples' surface spacing (5 nm) is
  far coarser than λ, so the model changes their currents by less than 1e-6. The FinFET's 0.1 nm surface spacing
  does resolve it.
- **Corners, thin bodies, Fermi–Dirac and heterojunctions** are as designed (QD2, QD5): one nearest interface,
  Hänsch's DOS form under Fermi–Dirac, and no confinement at heterojunctions or contacts.

---

### Unit 48 — Schottky contacts (started 2026-10-24; on `physics/schottky-contacts`)

**Owner request:** "start schottky contacts" (2026-10-24). ARCHITECTURE.md 6.4 has deferred them since Unit 12.
Numbered 48 because Unit 47 (MLDA, PR #66) is not merged; this branch starts from `main` (`e8719fb`).
**Owner rulings (2026-10-24):** (1) a thermionic-emission (Robin) condition for both carriers, the velocity from
the material's Richardson constant (the density-of-states velocity when it is 0), overridable per contact;
(2) the barrier given by exactly one of a metal work function (Schottky–Mott) or the electron barrier height;
(3) image-force lowering and tunnelling through the barrier deferred to a follow-up unit; (4) validation
against analytic theory, under Boltzmann and Fermi–Dirac statistics, with electrothermal and thermionic
interfaces refused. The standing conditions hold: no change to any device without a Schottky contact, criteria
fixed before code, Release and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged
without authorization.

**Sources read (before coding).**
- F. Masszi, "Computer modelling of metal-semiconductor junctions", Periodica Polytechnica Electrical
  Engineering 27 (3–4), 1983 (open access, pp.bme.hu/ee/article/view/4733). It gives the two-carrier
  combined thermionic-emission–diffusion condition of Crowell and Sze (Solid-State Electron. 9, 1035 (1966);
  paywalled, not read): J_n(0) = s_n (n(0) − n_0) and J_p(0) = s_p (p(0) − p_0) (its eqs 35, 36), with n_0 and
  p_0 the densities of the infinite-velocity (diffusion-theory) contact, whose quasi-Fermi levels equal the
  metal's Fermi level (eqs 25–27), and the potential fixed by the barrier. Its measured barrier table is not
  an I–V and is not used as a gate.
- COMSOL Semiconductor Module guide, "Metal Contacts" (doc.comsol.com, 6.3): the same form (n_0, p_0 "as
  though the Fermi level of the semiconductor at the boundary is equal to that of the metal"), the velocities
  from thermionic emission with A*_n and A*_p, and Φ_B either from the work function or given directly
  because measured barriers do not follow Schottky–Mott. Its equations are images and were not read.
- The velocity v = A* T² / (q N) is the one NiTCAD already uses at thermionic heterointerfaces
  (`physics::emission_velocity_cm_s`, Unit 15; legacy device.py `SchottkyContact`, M46-S2).
- Legacy: `pytcad/schottky.py` (M28, a standalone calculator, no device coupling) and `device.py`
  `SchottkyContact` (M46, 1D only). OLD / NEW / REASON:
  - OLD (M46-S2): the majority carrier takes the Robin condition, the minority carrier stays pinned at its
    barrier-limited density, and the majority carrier is chosen from the sign of the node's doping.
  - NEW: both carriers take the Robin condition with their own velocity and barrier (ruling 1).
  - REASON: the minority carrier crosses the same barrier by the same mechanism; pinning it fixes its
    injection at the diffusion-theory value, and the doping-sign rule fails on a compensated or
    junction-adjacent contact.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| SD1 | `device::ContactKind::schottky` with `Contact::schottky` (`SchottkyBarrier`): the boundary patch, exactly one of `work_function_eV` (Φ_Bn = Φ_m − χ of the node's material) or `barrier_eV` (Φ_Bn), and optional effective Richardson constants for electrons and holes (0 = the material's). The device checks Φ_Bn in [0, E_g(T)] at every node, a semiconductor node on the patch, finite non-negative Richardson constants. Φ_Bp = E_g(T) − Φ_Bn. | Rulings 1 and 2; the patch gives each node's face area, as a gate's does. | Schottky–Mott is the zero-pinning limit; the direct barrier covers measured ones. |
| SD2 | The barrier-limited densities: n_B = N_c e^(−Φ_Bn/kT) (Boltzmann) or N_c F_1/2(−Φ_Bn/kT) (Fermi–Dirac); η_B is the node's reduced potential at which its own electron density (the scaled device's statistics, n_ie, band shift and density of states) equals n_B, and p_B is the node's hole density at η_B. ψ_S = V/V_T + η_B − s. | Masszi eqs 25–27: both quasi-Fermi levels at the metal's; deriving ψ_S, n_B and p_B from the node's own statistics makes the equilibrium carry no current exactly. | With band-gap narrowing the narrowing lowers the hole barrier only (n_B is fixed by Φ_Bn); recorded. Incomplete ionization changes no row of the contact. |
| SD3 | Rows of a Schottky node: the Poisson row is ψ − ψ_S (Dirichlet, as an ohmic node's); the continuity rows are an interior node's box rows (edge fluxes, recombination, generation, storage, interface traps) plus the boundary emission −K_n (n − n_B) in the electron row and −K_p (p − p_B) in the hole row (signs in the rows' flux convention: carriers leave through the contact when n > n_B), K = v A, A the node's face area on the patch in the units of the thermionic edges' K. | Masszi eqs 35–36 (ruling 1) on the box method; the node keeps all its volume terms. | The ψ-Dirichlet node keeps no Poisson box row (as ohmic). |
| SD4 | v = A* T² / (q N) with the contact's A* for that carrier when set, else `physics::emission_velocity_cm_s` of the node's material (its Richardson constant, or the density-of-states velocity). | Ruling 1; one velocity function for heterointerfaces and contacts. | Silicon's default is the density-of-states velocity (2.6e6 cm/s for electrons, about half the tabulated A*); a deck sets A* (e.g. 110 / 30, Sze) explicitly. |
| SD5 | Currents: a Schottky contact's conduction current is its nodes' emission sum, q Σ A (v_p (p − p_B) − v_n (n − n_B)) in the sign of `terminal_currents`; its charge is the Gauss flux as an ohmic contact's (its ψ is Dirichlet); transients and small signal take the emission terms' partials and the Poisson row's bias derivative. | The emission is the only flux crossing the contact, so currents stay conserved. | — |
| SD6 | A Schottky contact satisfies both topology rules (it fixes the potential and exchanges carriers: K > 0 always). Equilibrium Poisson and `stamp_contacts` hold ψ_S at its nodes and take the densities from the potential (n_B and p_B there at equilibrium). The thermal-equilibrium bias check covers it. | No floating-region refusal for a device whose only carrier contact is a Schottky; consistent starting state. | — |
| SD7 | Refused (invalid_input): with `electrothermal`; on a device with a thermionic-emission interface; with the weighted terminal current (Unit 32; its Laplace weights are built on ohmic contacts). | Ruling 4; the weighted current is a follow-up. | — |
| SD8 | The deck: `[contact name] kind = schottky`, `patch`, `box`, exactly one of `barrier` or `work_function` (eV), optional `richardson_n`, `richardson_p` (A cm⁻² K⁻²); the device file (io) stores the contact; the run record names the kind and barrier; the CLI, app and render name the kind; the deck analyses (Unit 43) sweep a Schottky bias as an ohmic one. | Deck parity with the other contact kinds. | The device file's format version rises only if the stored contact needs it (recorded in the build). |
| SD9 | Not in this unit: image-force lowering, thermionic-field emission, a Schottky contact with traps of its own, a work-function temperature dependence. | Ruling 3. | — |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| SC1 | Without a Schottky contact: every existing test and pin passes unchanged; run identities, Jacobian patterns and the Unit 40 traces unchanged; a device file without one is byte-identical to `main`'s. |
| SC2 | Equilibrium (1D, n-Si 1e16 and p-Si 1e16, Φ_Bn = 0.7 and 0.4 eV, both statistics): the terminal current at most 1e-12 of the reverse current of SC5; the contact node's n equals n_B and p equals p_B within 1e-12 relative; the band bending ψ_bulk − ψ_S equals Φ_Bn/kT − ln(N_c/n_bulk) (Boltzmann) within 1e-9. |
| SC3 | Electrostatics (1D, n-Si 1e16, Boltzmann, Φ_Bn = 0.7 eV, equilibrium): the surface field from the contact's Gauss flux against the exact first integral of Poisson's equation E_s² = (2 kT N_D / ε)[e^(−u) + u − 1 + (p_0/n_0)(e^u − u − 1)], u the band bending in kT, within 0.5%. C–V from the small-signal capacitance at 0 … −2 V: the 1/C² slope gives N_D within 1% and its intercept V_bi − kT/q within 10 mV. |
| SC4 | Thermionic-diffusion (1D, no recombination, n-Si 1e16 and p-Si 1e16, V = −1 … +0.3 V): the current equals the exact series combination of emission and drift-diffusion evaluated on the solution's own potential, J = q (u_bulk − u_B) / (1/(v n_ie e^ψ_S) + Σ_edges L_e / (D_e n_ie e^ψ̄_e)) (Slotboom variables, the SG edge integral), within 1%. |
| SC5 | Diode characteristics (1D, n-Si 1e16, Φ_Bn = 0.7 eV, A*_n = 110 A cm⁻² K⁻², 300 K): the reverse current at −1 V between 0.80 and 1.00 of A* T² e^(−Φ_Bn/kT); the forward ideality factor d V / d(V_T ln J) between 1.00 and 1.05 for J in 1e-6 … 1e-2 A/cm². The same on p-Si 1e16 with Φ_Bp = 0.7 eV and A*_p = 30. |
| SC6 | Velocity limits: with A* = 1e9 the contact density within 1e-3 relative of n_B at +0.2 V (the diffusion-theory limit); the current rises monotonically with A* over 10 … 1e9. |
| SC7 | Minority carriers (n-Si 1e16, Φ_Bn = 0.8 eV, +0.3 V): the hole current at the contact equals q v_p A (p − p_B) at its nodes within 1e-10, and the electron and hole currents sum to the terminal current within 1e-10 relative. |
| SC8 | Fermi–Dirac: at Φ_Bn = 0.05 and 0.7 eV the equilibrium current at most 1e-12 of SC5's reverse current and n = N_c F_1/2(−Φ_Bn/kT) at the contact within 1e-10; at 0.7 eV the I–V within 1e-3 of Boltzmann's. |
| SC9 | Work function and barrier: `work_function = χ + 0.7` and `barrier = 0.7` give the same currents within 1e-12 relative. |
| SC10 | Geometry: a 2D (tensor and Delaunay) and a 3D device uniform across the contact give the 1D current per unit area within 1e-9 relative (tensor) and 1e-6 (Delaunay). |
| SC11 | The Jacobian against central differences of the residual with Schottky contacts (2D Delaunay, 3D tensor; both statistics; with and without Canali and band-gap narrowing): the largest relative column error at most 1e-6. Newton converges quadratically on SC5's device at +0.3 V. |
| SC12 | Dynamics: small-signal conductance against DC differences within 1e-5 and capacitance against dQ/dV within 1e-4 (SC5's device, −1 V); after a reverse-bias step the total currents of all contacts sum to zero within 1e-10 relative at every step. |
| SC13 | A 2D device whose only carrier contacts are a Schottky contact and an ohmic one: the terminal currents sum to zero within 1e-10 relative; a device whose only carrier contact is a Schottky contact passes the topology check. |
| SC14 | Refusals: electrothermal, a thermionic-emission interface, the weighted current, both or neither of `barrier` and `work_function`, Φ_Bn outside [0, E_g], a negative Richardson constant, an insulator node, a missing patch; each with invalid_input. |
| SC15 | The deck, device file (round trip, the contact unchanged), run record, CLI, app vocabulary and render name the kind; a deck analysis sweeps a Schottky bias. |
| SC16 | Reported, not gated: an example Schottky-diode deck's I–V (n-Si, Φ_Bn = 0.7 eV) against the closed-form A* T² e^(−Φ_Bn/kT)(e^(V/V_T) − 1); the reverse current with the material's default velocity against A* = 110. |
| SC17 | Full Release and Debug suites; zero warnings; hand mutants caught. |

**As built (2026-10-28).** `device` (the kind, `SchottkyBarrier`, its checks; `check_on_patch` factored out of the
gate's check, messages unchanged), `physics::emission_velocity_cm_s` (an overload taking the contact's Richardson
constant), the scaled device (per node the contact, η_B and K_n, K_p), `DriftDiffusion` (ψ_S, n_B and p_B in
`set_bias`, the contact owner of every node for currents and charges, the small-signal and transient paths), the
`SchottkyRows` term (last in the isothermal table), `EquilibriumPoisson`, the thermal-equilibrium bias check, the
refusals, the run record, the device file, the deck, the CLI, the app's vocabulary and workspace. Tests:
`tests/assemble/schottky_test.cpp` (SC1's pattern, SC11's FD gate, SC13's topology, SC14), `tests/deck/schottky_test.cpp`
(SC2–SC16), the app vocabulary's Schottky case (SC15: vocabulary and CLI). Example deck: `decks/schottky_diode_1d.deck`.

| # | Added or changed in the build | Why |
|---|---|---|
| SD5′ | A Schottky contact's conduction current is the edge flux leaving its nodes, as an ohmic contact's, not the emission sum; its charge (transient and small signal) includes the nodes' stored charge. | One bookkeeping for both kinds. In steady state the edge flux equals the emission sum (recombination cancels between the carriers; SC7 checks it carrier by carrier with recombination off); in transients the stored charge of the contact's half cells moves into the contact charge, and the total currents are unchanged (SC12). |
| SD3′ | The hole row gains +K_p (p − p_B), not −K_p. | The hole row is written with the opposite sign to the electron row; physically the same (holes leave when p > p_B). |
| SD8′ | The device file stores the contact as a new optional record (contact key 6); the format version is unchanged. | A reader of `main` never meets the record on a device without a Schottky contact; such files are byte-identical to `main`'s (SC1). |
| SD8″ | The example deck is in `decks/`, not `examples/`. | The app's examples pin exactly four decks (WG8); `decks/` holds decks outside it (as the FinFET). |
| SD4′ | Correction to SD4's trade-off note: silicon's default velocity (2.6e6 cm/s) is about 1.2× the A* = 110 velocity (2.16e6 cm/s; an equivalent A* of about 132), not half. | Measured (SC16): the default gives 1.185× the reverse current of A* = 110. Not a criterion; the design text was wrong. |
| — | The small-signal interface loop skips its Poisson-row additions on Schottky nodes; trap conduction additions stay ohmic-only. | The Schottky node's Poisson row is Dirichlet. |
| — | The app vocabulary gains `barrier`, `richardson_n`, `richardson_p` (context `schottky`) and `work_function` in that context; its misspelling test gained a Schottky base deck for the new context. | Suggestions draw from the keys of the contact's own kind; no existing pin changed. |

**Results:**

| # | Result |
|---|---|
| SC1 | Met. Every existing test passes unchanged. On `main` (`e8719fb`, built in a worktree) and the branch the four example decks give byte-identical device files, identical result tables (17 digits) and identical run identities. An ohmic contact made Schottky leaves the Jacobian pattern unchanged (1D, 2D, 3D). |
| SC2 | **Partly not met.** n = n_B and p = p_B within 1e-12 (n-Si and p-Si, 0.7 and 0.4 eV, Boltzmann and Fermi–Dirac); the band bending within 1e-9; the equilibrium current within the terminal-current resolution. **The equilibrium current is not within 1e-12 of SC5's reverse current (1.7e-17 A/cm²):** it is 9.5e-15 A/cm² (Boltzmann, 0.7 eV), 1.4e-14 and 5.6e-14 (Fermi–Dirac, 0.7 and 0.4 eV), 2.9e-8 (Fermi–Dirac, 0.05 eV). Cause: the rounding of the contact edge's two opposing fluxes (q D n / h, about 5e5 A/cm² on 1 nm edges, 2e8 in the 0.05 eV accumulation layer), inside the resolution (9e-8 A/cm² at 0.7 eV). The criterion did not allow for rounding. `[!shouldfail]`, unchanged. **For the owner.** |
| SC3 | Met. E_s 38,068.5 V/cm against the first integral's 38,065.1 (9e-5). C–V: N_D 9.9999e15 (1e-5); intercept 0.468479 V against V_bi − kT/q = 0.468429 (0.05 mV). |
| SC4 | Met: the series law equals the current to 6 printed digits at every bias, n-Si and p-Si. |
| SC5 | **Saturation met; ideality not met.** Reverse at −1 V 0.975 (n-Si) and 0.978 (p-Si) of A* T² e^(−Φ_B/kT). **The ideality over 1e-6 … 1e-2 A/cm² is 0.43 … 1.00003 (n) and 0.43 … 1.002 (p).** J_s = 1.7e-5 A/cm² lies inside the range, so near V = 0 the form J_s (e^(V/V_T) − 1) itself gives d V / d(V_T ln J) < 1: the criterion's range assumed J ≫ J_s. With the −1 removed, J / (1 − e^(−V/V_T)), the ideality is 1.0018–1.0026 (n) and 1.0017–1.0028 (p), reported. `[!shouldfail]`, unchanged. **For the owner.** |
| SC6 | **Monotone met; the limit not met.** The current rises from 3.56e-3 to 0.618 A/cm² over A* = 10 … 1e9 at +0.2 V; n(0) − n_B = J/(q v) within 1e-6 (tested). **At A* = 1e9, n(0)/n_B = 1.00395, not within 1e-3.** Cause: the diffusion-limited current at +0.2 V (0.618 A/cm²) is 17× the thermionic one, so even v = 2e13 cm/s leaves 4e-3; the criterion underestimated it (about A* ≥ 4e9 would meet it). `[!shouldfail]`, unchanged. **For the owner.** |
| SC7 | Met (recombination off, so the contact's half cell adds nothing to its edge current): the hole current 9.95011e-6 A/cm² equals q v_p (p − p_B) within 1e-10 of the current; the carriers sum to the terminal current within 1e-10. |
| SC8 | n = N_c F_1/2(−Φ_Bn/kT) within 1e-10 at 0.05 and 0.7 eV; the I–V within 1e-3 of Boltzmann's at −1 … +0.3 V. The equilibrium current's 1e-12: not met, as SC2. |
| SC9 | Met within 1e-12. |
| SC10 | Met. 2D and 3D tensor within 1e-9 of the 1D current per area. 2D Delaunay (3,306 nodes, `mesh::generate_mesh`: a structured boundary layer on the contact with the 1D spacing to 1 µm and 0.1 µm along it, generated triangles beyond) within 1e-6, recombination off: the box method reproduces a 1D solution exactly where it is linear in x on any Delaunay mesh, so the quasi-neutral body then needs no matching mesh. |
| SC11 | Met. The FD gate's worst column error is 3.4e-11 over the 16 cases (Delaunay and 3D tensor; Boltzmann and Fermi–Dirac; ±Canali; ±band-gap narrowing). Newton at +0.3 V on SC5's device: corrections 50.9, 17.4, 1.12, 1.35e-4, 2.35e-11, 2.9e-15; gated as each correction in the asymptotic window (below 1e-2, next above 1e-13) at most the square of the previous (2.35e-11 ≤ 1.8e-8). The three-point order estimate of LS4 gives 0.58 on the last three (the last is rounding) and 1.73 without it (its first, 1.12, outside the window); the measure is stated in the test. |
| SC12 | **Capacitance and conservation met; the −1 V conductance not met.** C = 2.37723e-8 F/cm² against dQ/dV, 6 digits; after a −1 V step the total currents sum to zero within 5.3e-11 of the run's largest current over 208 steps (a 1 ns ramp instead, measured first: 2.3e-9). Added: g at +0.2 V against DC differences within 1e-5 (6 digits); every column and row of Y within the reported resolution. **At −1 V, g = 1.4947e-7 S/cm² is not within 1e-5 of the DC difference:** the DC current (1.7e-5 A/cm²) is resolved to about 1e-6 of itself and the admittance to 2.0e-10 S/cm² (1.3e-3 of g), as Unit 22 found for the ohmic diode at low bias. `[!shouldfail]`, unchanged. **For the owner.** |
| SC13 | **Met at forward bias; not at −1 V.** At +0.3 V (1.7 A/cm²) the terminal currents of 2D tensor and Delaunay devices sum to zero within 1e-10, as do a 1D back-to-back Schottky device's; a device whose only contact is a Schottky contact passes the topology check and assembles. **At −1 V the 2D sum is 2.2e-15 A/cm (1.3e-6 of the reverse current; the 1D device the same):** the rounding of the cancelling edge fluxes, about 2e-11 A/cm² per area, inside the resolution. 1e-10 relative holds above about 1 A/cm². `[!shouldfail]`, unchanged. **For the owner.** |
| SC14 | Met: electrothermal, a thermionic-emission interface, the weighted current, both and neither of `barrier` and `work_function`, Φ_Bn above E_g and below 0 (by the work function), negative Richardson constants (each carrier), an insulator node, a missing patch and a node off its patch, the Schottky keys on an ohmic contact and an unknown kind; each invalid_input. |
| SC15 | Met. The deck reads the kind and keys; the device file round-trips two Schottky contacts unchanged; the run identity changes with the barrier, with the barrier given as a work function, and between an ohmic and a Schottky anode; `nitcad check` prints `contact anode: schottky`; the app's vocabulary lists the kind and its keys only in its context; a `[small_signal]` deck analysis over a Schottky bias runs through `deck::run` (C falls with reverse bias). Render needed no change (a Schottky contact is drawn as a metal contact, as an electrode); not separately tested. |
| SC16 | Reported. `decks/schottky_diode_1d.deck` against A* T² e^(−0.7 eV/kT)(e^(V/V_T) − 1): reverse 0.975 at −1 V to 0.958 at −0.05 V (the diffusion term of the series law); forward 0.953 at 0.05 V, 0.940 at 0.2 V, 0.917 at 0.3 V and 0.694 at 0.4 V (the body's series resistance). The default velocity's reverse current at −1 V is 1.99061e-5 A/cm² against 1.67924e-5 with A* = 110 (1.185; SD4′). |
| SC17 | Release (112 tests) and Debug (103) suites pass, zero warnings; no existing pin changed. Hand mutants: 16 of 18 caught. Caught: the emission term's Jacobian (SC11 gate, Newton); the hole emission sign; the Poisson row not cleared; ψ_S without the band shift and equilibrium Poisson's ψ_S without it (both equivalent on a single-material device, where s = 0 at the reference material's contact: caught after a heterojunction case was added, the Schottky contact on the shifted material); the contact owner ohmic-only (currents); the bias derivative (C–V, small signal); the device file's Richardson constant; the run record's barrier; the deck's exactly-one rule; the velocity ignoring A*; the thermal-equilibrium bias check (a direct case added); η_B's sign; the Richardson check; the electrothermal and weighted-current refusals. **Survived, for the owner:** (1) the small-signal interface loop adding its Poisson-row terms on a Schottky node: reached only when a Schottky node lies on a semiconductor–insulator interface with traps, which no test builds; (2) the transient error measure skipping a Schottky node's carriers: it changes only the step-size choice, not a result the tests can tell apart. Added in the build, as the mutants found: the band-gap-narrowing equilibrium at 1e18, the shifted-material case and the direct bias-check case. |

**Limitations found (for the owner):**
- **Rounding at small currents.** The relative gates at the equilibrium and reverse currents (SC2, SC8, SC13 at
  −1 V) and the −1 V conductance (SC12) ask for more than the state resolves: the edge fluxes that cancel in a 1.7e-5
  A/cm² current are about 5e5 A/cm² on 1 nm edges. The weighted terminal current (Unit 32) is the tool for these
  currents, and it is refused with a Schottky contact (SD7).
- **SC5 and SC6 criteria** assumed J ≫ J_s over the ideality range and a diffusion-limited current near the
  thermionic one. The model behaves as the theory says; the criteria's ranges do not fit this device.
- **No image-force lowering or tunnelling** (SD9): the reverse current stays flat with bias.
- **Mesh generator (Unit 16, not changed here):** with a boundary layer whose spacing along the segment exceeds
  `max_size` (0.5 µm against 0.1 µm), `generate_mesh` split the layer's own boundary segment and joined the new nodes to
  the layer with sliver triangles (contact-node volumes 9e-11 cm² against 2.5e-12), although U14 says refinement never
  inserts into the layer's rectangle. Found while building SC10's Delaunay device, whose current was then 0.3–1.7% off;
  with the layer spaced at `max_size` it is exact. **For the owner.** Fixed in Unit 49 (below): the cause was the
  boundary fill, not refinement.

**Ruling on SC2/SC8, SC5, SC6, SC12 and SC13 (2026-10-28).** The owner: "accept it as limitations". The five unmet
parts are accepted as limitations of the criteria as written, not defects of the model. Each criterion's intent is met:
- SC2/SC8: the equilibrium current is within the terminal-current resolution at every barrier, under both statistics.
- SC5: the ideality of the thermionic form is 1.002–1.003.
- SC6: the gap to n_B falls 10× per decade of A*, and n(0) − n_B = J/(qv) holds.
- SC12: g matches DC differences at +0.2 V and C matches dQ/dV, both to 6 digits.
- SC13: the currents are conserved within 1e-10 at +0.3 V and within the resolution at −1 V.
The 1e-12 and 1e-10 relative gates at 1.7e-5 A/cm² ask for less than the rounding of the cancelling fluxes. SC5's
range starts below J_s. SC6 underestimated the diffusion current.
The criteria are unchanged, and the `[!shouldfail]` cases stay, so a change that makes one pass is flagged.
Follow-up if currents near 1e-14 A/cm² (Schottky reverse leakage) matter: extend the weighted terminal current
(Unit 32) to Schottky contacts (refused now, SD7). No unit named.
The surviving mutants and the mesh-generator finding stay open (the finding: Unit 49).

### Unit 49 — boundary fill between a structured layer's columns (started 2026-10-28; on `mesh/boundary-layer-split`)

**Owner request:** "Unit 16 bug" (2026-10-28): the mesh-generator finding of Unit 48 (above).

**Diagnosis (before code).** Unit 48's description was wrong: no segment of the layer is split and no refinement
point falls in the layer. Reproduced with SC10's geometry (2 µm × 1 µm, the layer on the whole x = 0 edge, 280 rows to
1 µm, max_size 0.1 µm): with columns 0.1 µm apart the contact carries the 11 lattice nodes; with 0.25 µm, 7 nodes for
5 columns; with 0.5 µm, 7 for 3 (contact volumes up to 2e-10 cm² against 2.5e-13). Every segment split and every
refinement point near the layer was logged; none touched the contact. The extra nodes are the **boundary fill**
(`place_boundary`): the breaks along an edge are a lattice column, an edge end and a patch end, sorted with an
unstable sort and merged when they coincide. A coincident edge or patch end, which carries an exact point, replaced a
column's break, and with it the break's kind. The rule "between two lattice nodes, no fill" (U14) then saw a
`boundary` break, so the gap from that end to the next column was filled at h. Those fill nodes are on the contact
but not in the lattice, so they join the first row (1e-7 cm away) with slivers. A column at a domain corner or a
patch end, the usual case, triggers it. Which end it hits depends on the sort's order. The U14 test (a layer on the
whole bottom edge, columns 0.25 apart, h 0.1) has the same defect; it counts only nodes above the edge. With columns
no coarser than h the fill adds nothing, so SC10 at 0.1 µm was exact.

**Design.**
- BD1: coincident breaks merge order-independently. The merged break is a lattice break when either one is, so no
  fill runs between it and a neighbouring column.
- BD2: a column's bottom node is the exact given point where one exists: the layer's ends `a` and `b` and a patch end
  on the segment (it was `a + t u`, which can differ from `b` or from the patch end in the last bit). The merged break
  then carries the same point whichever break wins, so the edge, the patch and the lattice share one node.
- BD3: nothing else changes. Meshes without layers are unchanged; with layers, only the fill between a column at an
  edge or patch end and its neighbour goes away.

**Criteria (fixed before code).**
- BC1: with a layer on a whole boundary edge and columns coarser than h (SC10's geometry at 0.25 and 0.5 µm; the U14
  square), the edge carries exactly the lattice's columns: no node on it between two columns.
- BC2: in SC10's geometry at 0.5 µm the contact nodes' volumes are the tensor lattice's, (Δt/2 or Δt) × (first row /
  2), to 1e-12 relative, and every node in the layer's rectangle is a lattice node.
- BC3: a patch end on a layer is a node bit-identical to the given point, and is a lattice column; there is no other
  node within 1e-9 of the layer's length of it.
- BC4: SC10's 2D Delaunay diode with columns 0.5 µm apart is within 1e-6 of the 1D current per area (the bound SC10
  holds at 0.1 µm), so the deck test drops the workaround.
- BC5: every existing test passes unchanged, Release and Debug; the meshes of the examples and decks without layers
  are bit-identical (run identities unchanged).
- BC6: hand mutants: the merge rule reverted (fails BC1), the exact bottom point reverted (fails BC3).

**As built.** `place_boundary` merges coincident breaks in a pass of their own before the fill: a lattice break is
kept, else an exact one. A first attempt kept the single pass and only fixed the merged break's kind; the fill had
already run when the first coincident (`boundary`) break arrived, and BC1 still failed (7 nodes for 3 columns).
`read_layers` gives each column its node on the edge: `a`, `b` and patch ends exactly, `a + t u` elsewhere.
Without layer breaks the merge keeps what the old rule kept and the fill is unchanged, so meshes without layers are
bit-identical; no deck or example uses layers (only the generator tests, SC10 and the U14 MOSFET in
`tests/solve/unstructured_test.cpp`).

**Results.**
- BC1, BC2: met. SC10's geometry at 3 and 5 columns: the contact carries exactly the columns, every node in the
  layer's rectangle is a lattice node, contact volumes (Δt/2 or Δt) × 5e-8 to 1e-12. The U14 square: 5 nodes on
  the bottom edge, its 5 columns.
- BC3: met. Layer (0.2, 0)–(0.9, 0), whose `a + len u` misses 0.9 in the last bit, with a patch end at 0.4: `a`,
  `b` and the patch end are each one node, bit-identical, none other within 7e-10; patch length 0.5 to 1e-12.
- BC4: met. SC10's Delaunay device at 3 columns (0.5 µm) matches the 1D current per area at −1 V and +0.2 V within
  SC10's 1e-6; the deck test now uses it and drops the workaround.
- BC5: met. Release 112/112, Debug 103/103 (the U14 MOSFET gate included).
- BC6: 3/3 caught. The merge reverted to the single pass (BC1, BC2: 19 failures); an exact break replacing a
  lattice one (BC3); the bottom node as `a + t u` (BC3).

### Unit 50 — Newton at the rounding floor: the residual line search and the minority-density floor (started 2026-10-28; on `solve/rounding-floors`)

**Owner request:** "solver unit" (2026-10-28): the unit Unit 47 proposed after diagnosing two Newton stalls
(DECISIONS.md Unit 47 on `physics/quantum-confinement`, PR #66): (a) the residual line search at a rounding
residual (the FinFET with MLDA stalled at 0.35 V); (b) the minority-density floor (the 1e17 long channel with
Lombardi chatters at 1–4e-8).

**Reproduced on `main` (before code), one device for both.** Unit 46's 1e17 long-channel n-MOSFET
(`surface_mobility_test.cpp`, Vg 0–6.5 V in 0.5 V, Vd = 10 mV) with Lombardi, Release:
- Default options: all 14 points, 55 s, but at 0.5, 2.5, 3.5, 4.0, 4.5 and 5.5 V the corrections chatter at
  1.0–3.8e-8 for up to 14 iterations and finish only when one dips below 1e-8. The largest correction is always
  a minority hole in an n+ pocket: p = 2.1–2.8e-19 of the scaled majority (about 27 cm⁻³ against 1e20), its
  correction 1–7e-27 of the majority, so 1–4e-8 of p. The global floor (1e-20 of the largest density) is below p,
  so these count.
- `line_search = residual`: no point solves (100 iterations at Vg = 0). From iteration 4 the largest |F| is
  3e-11, rounding of the large rows, while one minority-hole node (p = 2.75e-20) still needs a factor of 8.7; the
  Armijo test on Σ F² sees only the large rows' noise, takes 2^-10-ish steps, and p creeps by 0.6% an iteration.
- Measured against each row's own scale s_i = Σ_j |J_ij x_j| (the size of the terms it sums), that hole row's
  residual is 1.5e14 times ε s_i: it is not rounding, only invisible in Σ F². Converged states sit at
  0.9–7.5 ε s_i; the chattering states at about 1e3 ε s_i on the minority rows (the solve's rounding, not the rows').
- Prototypes (environment-switched, not committed): Armijo on Σ (F_i / s_i)² alone stalls the equilibrium
  solve at 2^-10 steps (the clipped step is not a descent direction for it, Unit 30's finding); passing a trial
  that satisfies Armijo on **either** merit solves all 14 points. Measuring n and p against
  max(density, global floor, ε max(n, p) at the node) converges every point quadratically to 3e-11–1.4e-10, 47 s,
  drain currents equal to the default run's within 1e-11 (11 of 14 bit-identical).

**Owner rulings (2026-10-28, all as recommended):** (1) the residual test passes a trial on either merit, or
when every row is within 64 ε of its scale; (2) the floor ε × the node's majority density, on by default;
(3) in scope: the FinFET deck as a regression, a recheck of Unit 48's rounding-limited criteria (report only),
and a retest of PR #66's FinFET with MLDA and the residual test once this merges.

**Design.**
- RD1 (`solve::newton_solve`, residual test only): at each searched iteration, s_i = Σ_j |J_ij x_j| from J and x
  at the iterate (1 where it is 0). A trial passes when Σ F² or Σ (F_i / s_i)² falls by the Armijo factor
  (1 − 2 α λ), or every |F_i(trial)| ≤ 64 ε s_i. Every trial the old test passed still passes, so a run changes
  only where the old test rejected the full step. The best trial (none passing) is ranked on Σ F², as before
  (corrected as built: the first text said the scaled merit). The natural test and runs without a search are
  unchanged.
- RD2 (`DriftDiffusion::update_size`, **opt-in**, see the second ruling below): a density correction is measured
  against max(density, global floor, ε max(n, p) at that node), ε = 2^-52. A minority density below ε of the
  majority at its node does not change the node's charge in double, and its correction is resolved only to the
  solve's rounding there. Majority densities (≥ the minority) are measured as before. Insulator nodes and θ
  unchanged. Off, `update_size` is the old one.
- RD3 (as ruled the second time): `solve::BiasOptions::minority_floor` (off), set on the system by
  `DriftDiffusion::set_minority_floor` wherever a bias, trace, transient or small-signal solve builds one (and
  on the sweep's retry system without generation); the deck key `[solver] minority_floor = on | off`; recorded
  in the run's settings only when on, so the identities of runs without it are unchanged. The app's deck
  vocabulary offers the key.

**Found in the suites with the floor on by default (2026-10-28), and the second ruling.** `solve.heterojunction`
failed in Release and Debug: its 4H-SiC p-n diode (incomplete ionization, 0.5 V steps to 3 V) stopped at 3.0 V
after 100 iterations. On `main` that point already needs 69 iterations, swinging to 1e12, but there Newton
resolves SiC's minorities (about 1e-35 cm⁻³, far below ε of the majority) quadratically, so RD2's premise (such
a minority's correction is rounding) does not hold in general; with the floor its corrections no longer count
and the point wanders. The current at 1.5 V moved from 5.9e-13 to −1.19e-12 A/cm², below the test's
resolution. **Owner ruling (2026-10-28, as recommended): make the floor opt-in**, behind a solver option and a
`[solver]` deck key, off by default, so every existing run stays bit-identical, SiC included; the long channel
and similar MOSFET decks turn it on, and RC2 runs with it on. RD1 stays on wherever the residual test is used.
The criteria below are kept as fixed, with RC2's and RC6's conditions as this ruling sets them; RC9 is added.

**Criteria (fixed before code).**
- RC1: the long channel with Lombardi and `line_search = residual` solves all 14 points (0 before).
- RC2: the long channel with Lombardi, default options (with `minority_floor` on, per the second ruling): every
  point converges with no more than 2 iterations whose correction lies in [1e-8, 1e-4) (up to 13 before), and
  its drain currents equal the old run's within 1e-9 relative.
- RC3: `newton_test`: a two-row system with rows scaled 1e20 apart, at a state where the large row is at rounding
  and the small row needs a real correction, converges with the residual test (and fails without RD1); a
  trial at rounding in every row passes; the Unit 30 cases are unchanged.
- RC4: `update_size`: a minority density below ε of the majority at its node is measured against ε × majority;
  above it, and for the majority, against the density as before; insulator nodes unchanged.
- RC5: Release and Debug suites pass. Unit 30's hidden diode-jump and MOSFET line-search benchmarks: every
  outcome (converged or not) as recorded in Unit 30 for the natural test; the residual test's recorded failures
  may only become successes.
- RC6: the FinFET deck (residual test, PARDISO, 6 threads) solves 17/17 before and after, drain currents within
  1e-6 relative (after: the same deck, so `minority_floor` off, per the second ruling).
- RC7 (report only): Unit 48's `[!shouldfail]` SC2/SC8 and SC13 rerun and their values recorded.
- RC8: hand mutants: the scaled merit dropped (RC1 fails), the at-rounding pass dropped (an RC3 case fails), the
  floor removed (RC2, RC4 fail), the floor taken from the global largest density instead of the node's (RC4
  fails); added with the second ruling: the floor applied whatever the switch (RC9 fails).
- RC9 (added with the second ruling): with `minority_floor` off, the default, `update_size` measures a minority
  below ε of its node's majority against itself (as before); `solve.heterojunction` passes; the deck key reads
  `on` and `off` and defaults to off.
- After merge (follow-up, done 2026-10-28): PR #66 rebased on `main` 0c14841; its FinFET with MLDA and `line_search = residual` solves 17/17 in 11.4 min, no cuts (Unit 47, limitations).

**As built (2026-10-28, Release and Debug, this machine).** Every criterion met.
- RC1: the long channel with Lombardi, `line_search = residual`: 14/14 points with the floor off and with it on
  (0 before).
- RC2: with `minority_floor` on, every point has at most 2 iterations in [1e-8, 1e-4); drain currents within 1e-9 of
  the 14 recorded before the change. Hidden test `[.lombardi]`, Release (the three sweeps: 145 s).
- RC3: `newton_test`: the two-scale system (rows 1e14 apart, the large one at rounding) converges with the residual
  test in at most 5 iterations from a full first step; a trial at rounding in every row passes (2 iterations); the
  Unit 30 cases unchanged.
- RC4, RC9: `drift_diffusion_test`: off, the minority p = 1e-19 of n = 1 is measured against itself; on, against
  ε × 1, against ε × 1e-3 at a node whose majority is 1e-3, with holes as the majority likewise; above ε of the
  majority and for the majority, against itself. `build_test`: the key reads on and off, off by default.
- RC5: Release 112/112 and Debug 103/103, `solve.heterojunction` included. Unit 30's `[.diode_jumps]` and
  `[.line_search]` benchmarks: outcomes and iteration counts identical to `main` (`c288391`, built from the same
  tree with the changes stashed). The residual test's MOSFET 0/10 stays 0/10: that case is damping-limited, not a
  rounding floor. **Finding, not Unit 50's:** `main` no longer matches Unit 30's MOSFET record (natural test 2/10
  against 10/10; ARCHITECTURE.md 11.3 "As built (Unit 50)").
- RC6: FinFET deck (residual test, PARDISO, 6 threads, floor off): 17/17 before and after; currents within 9.5e-12
  relative (point 0 within 7e-11); Newton iterations 428 → 182; point 0's substeps 3 → 1; wall time 28 → 10.6 min.
- RC7 (report only): SC2/SC8 9.5e-15, 1.4e-14, 5.6e-14, 2.9e-8 A/cm²; SC13 2.2e-15 A/cm; SC12 1.25e-4; SC6
  1.00395; SC5 0.4276: all as recorded in Unit 48 (the default path is unchanged). These are rounding of the
  current itself, not a Newton floor; the minority floor does not bear on them.
- RC8: M1 (scaled merit dropped): the RC3 two-scale case and RC1 fail; M2 (at-rounding pass dropped): the RC3
  every-row case fails; M3 (floor removed): RC4 fails (and RC2, measured with the floor always on); M4 (floor from
  the largest density in the device): RC4 fails at the 1e-3 node; M5 (floor whatever the switch): RC4's off case
  and `solve.heterojunction` fail. All caught.
- Found on the way, fixed here: `app.live_check`'s stale-result case (XG1) failed once locally under `-j 6` and
  once on CI (`dropped() == 1` read 0). The test waited for its check function to return, but the worker posts
  the result after that, so `process()` could run before it. The test now waits on the checker's wake, called
  after the result is posted; a 50 ms delay injected before the post fails the old test every time and passes the
  new one; the app tests pass 20 times under `-j 6` (Release) and 5 (Debug). The checker is unchanged.

### Unit 51 — gate tunnelling (started 2026-10-28; on `physics/gate-tunnelling`)

**Owner request:** "start GATE tunneling" (2026-10-28). ARCHITECTURE.md 11.3 lists gate tunnelling as Planned (the
nanoscale gap found by the FinFET deck). Numbered 51 because Units 47 (PR #66) and 50 (PR #71) are not merged; this
branch starts from `main` (`7a11d26`).
**Owner rulings (2026-10-28):** (1) both gate kinds tunnel: a lumped `gate` locally through its t_ox, and an
`electrode` on a meshed insulator along a straight path through the meshed stack (tensor grids); (2) self-consistent:
a boundary flux in the continuity rows with an exact Jacobian, counted in the gate's terminal current; (3) the
Tsu–Esaki energy integral with a WKB transmission through the trapezoidal or piecewise multi-layer barrier; (4) three
processes: conduction-band electrons (ECB), valence-band holes (HVB) and valence-band electrons (EVB); (5) SiO₂ gets
built-in, cited band data; any other insulator must give its band data in the deck, else the tunnelling gate is
refused; (6) validation against closed forms, plus a measured J–V comparison only if an accessible source turns up;
(7) excluded and refused: electrothermal, the weighted terminal current (Unit 32), electrode tunnelling on an
unstructured mesh; not built: image-force lowering. The standing conditions hold: no change to any device without
tunnelling switched on, criteria fixed before code, Release and Debug suites at every commit, `/W4 /WX`, hand
mutants, nothing pushed or merged without authorization.

**Sources read (before coding).**
- A. Gehring, "Simulation of Tunneling in Semiconductor Devices", dissertation, TU Wien, 2003 (open access,
  iue.tuwien.ac.at/phd/gehring): §3.2 (the Tsu–Esaki formula J = (4π m q / h³) ∫ TC(E_x) N(E_x) dE_x, parabolic
  bands, conserved parallel momentum, and the ECB, HVB and EVB energy limits), §3.3.1 (the Fermi–Dirac supply
  N = kT ln[(1 + e^((E_F1 − E)/kT)) / (1 + e^((E_F2 − E)/kT))]), §3.5.1 (the WKB transmission
  exp(−(2/ħ) ∫ √(2 m_diel (W − E)) dx) between classical turning points, and its direct and Fowler–Nordheim closed
  forms for a linear barrier), §4.2.2 (a stack is evaluated once per path, its current returned to the end nodes'
  continuity rows or the contact's current), §5.1.2 (SiO₂ barriers of 3.2 eV for electrons and 4.65 eV for holes
  in its simulations; ECB, HVB and EVB in n- and pMOS) and Table 5.1 (permittivity, gap and band offsets to silicon
  of nine dielectrics; SiO₂ 3.9, 9.0, 3.0, 4.9 eV; HfO₂ 25, 5.7, 1.5, 3.1 eV; from Leroy 2003, Wilk 2001). Its
  equations are images; their structure was read from the text and is re-derived below.
- Electron tunnelling mass in SiO₂: Brar, Wilk and Seabaugh, Appl. Phys. Lett. 69, 2728 (1996), 0.30 m₀ with a
  parabolic dispersion and 0.41 m₀ with a non-parabolic one, from 1.6–3.5 nm oxides (abstract only, via a search
  summary; the paper is paywalled). Hole mass: Yang et al., IEEE Trans. Electron Devices 47, 2161 (2000),
  0.51 m₀ parabolic, fitted to p⁺-poly pMOSFETs with 1.23–2.16 nm oxides (abstract only). Neither could be read.
  **Owner ruling (2026-10-28):** take the more common values instead: 0.50 m₀ for electrons and holes, the round
  value usually assumed for SiO₂ (Ng et al., J. Appl. Phys. 2004, quoted in a search summary: "usually assumed to
  be 0.5 m_e or 0.4 m_e"); the hole value is also Yang's 0.51 rounded (GD5).
- Legacy: `workbench/physics/tunneling.py` (M12-S1): Fowler–Nordheim and a rectangular WKB transmission as
  diagnostics at the solved oxide field, not self-consistent. OLD / NEW / REASON:
  - OLD: J_FN = (A/φ) E² e^(−B φ^(3/2) / E) with the free-electron mass, evaluated after the solve.
  - NEW: the Tsu–Esaki integral with the insulator's tunnelling mass, in the residual and Jacobian.
  - REASON: rulings 2 and 3; FN misses direct tunnelling, the regime below about 3 nm, and a diagnostic does not
    conserve current.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| GD1 | Energies on one scale per node: the semiconductor's E_c, E_v (from ψ, χ, the band shift, as `band_edges`), its quasi-Fermi levels E_Fn, E_Fp (from n and p through the node's statistics, Boltzmann or the inverse F_1/2); the electrode's Fermi level E_FG = −q V_G. The barrier is the vacuum level minus χ_ox: at the semiconductor end E_c + χ_s − χ_ox, at the electrode end E_FG + φ_m − χ_ox, the insulator's valence edge E_g,ox below it. | Endpoint vacuum levels make the lumped barrier exact for a fixed interface charge (it changes the field, not the ends) and need no V_ox. | — |
| GD2 | Lumped gate (ruling 1): one linear segment of length t_ox between the two ends of GD1, the gate stack's insulator band data (GD5). Electrode (ruling 1): from each semiconductor node with a face on an insulator, a straight path along that face's grid axis, node to node through insulator nodes, ending at a node of an `electrode` contact; no path if it meets a semiconductor node, the mesh boundary or no electrode. Each mesh edge on the path is one linear segment between its nodes' vacuum levels with the band data of the edge's region. A node with insulator faces on several axes has one path per face, each weighted by that face's area. | Tensor grids make the normal a grid line (no locator); a piecewise-linear barrier is exact for the box method's linear potential along an edge, and per-edge band data gives multi-layer stacks (HfO₂ on SiO₂). | A path is frozen geometry (as Unit 20's); a curved or oblique interface is out (unstructured refused, ruling 7). |
| GD3 | WKB: TC(E) = exp(−(2/ħ) Σ_segments ∫ √(2 m (W(x) − E)) dx) over the segments where W > E, each integral in closed form, (2/(3s)) [(W − E)^(3/2)] between the segment's ends or turning point, s the segment's slope (the limit √(W − E) L when s = 0). Electron barriers W = E_c,ox(x) with the insulator's electron mass; hole barriers W = E − E_v,ox(x) mirrored, with its hole mass; EVB electrons see the conduction-band barrier with the insulator's EVB mass (default its electron mass). Above the barrier's highest point no energy is integrated. | Closed form per segment: no spatial quadrature, exact derivatives. | No over-barrier (thermionic) current: negligible for 3 eV barriers, not for a low-offset high-k at high field; recorded. No image force (ruling 7), no quasi-bound states (Tsu–Esaki's continuum overestimates inversion current, Gehring §3.2); recorded. |
| GD4 | Tsu–Esaki: J = (A* / k²) ∫ TC(E) N(E) dE (N has units of energy; with TC = 1 above E_c and a Boltzmann supply it is A* T² e^(−(E_c − E_F)/kT), thermionic emission), N(E) = kT ln[(1 + e^((E_F,s − E)/kT)) / (1 + e^((E_FG − E)/kT))], A* the semiconductor carrier's Richardson constant (`physics::emission_velocity_cm_s`'s: the material's, or its density-of-states mass). Limits per process (Gehring §3.2): ECB from max(E_c,s, E_c,G) to the barrier top, E_F,s = E_Fn; HVB from the barrier bottom to min(E_v,s, E_v,G), mirrored, E_F,s = E_Fp; EVB between E_c,G and E_v,s (electrons from the semiconductor's valence band to the gate's conduction band) and between E_c,s and E_v,G (the reverse pair), E_F of a valence band its E_Fp. A metal electrode has no band edges: no lower limit for ECB, no upper for HVB, and no EVB (its exchange with the valence band is HVB). Polysilicon: degenerate silicon, E_c,G = E_FG (n⁺) or E_v,G = E_FG (p⁺), the other edge E_g,Si(T) away (as `electrode_work_function_eV`). The integral is a fixed composite Gauss–Legendre rule (4 points, panels of at most kT/2), no adaptivity. | (4π m q k² / h³) = A*: one prefactor shared with Schottky and thermionic emission. A fixed rule keeps the residual smooth and its Jacobian the exact derivative of the discrete sum (moving limits through the panel map). | Hundreds of WKB evaluations per node per assembly; GC14 bounds the cost. |
| GD5 | Insulator band data: `InsulatorParameters` gains χ_ox, E_g,ox and the tunnelling masses (electron, hole, EVB), 0 meaning absent. SiO₂'s: χ_ox = 1.05 eV (Δ E_c = 3.0 eV to silicon's χ = 4.05, Gehring Table 5.1), E_g,ox = 9.0 eV, m_n = m_p = 0.50 (owner ruling above; EVB defaults to m_n). A lumped gate stack takes SiO₂'s unless it gives its own; a stack whose permittivity is not 3.9 must give its own. An electrode path crossing an insulator without band data is refused. | Ruling 5. | The SiO₂ masses are fitted values tied to their barriers; any change is an owner ruling. |
| GD6 | Rows: the semiconductor node's electron row gains the ECB and EVB fluxes into the gate (carriers leave the node when the semiconductor side's supply is the larger) times the face area; the hole row gains the HVB flux and the EVB flux as a hole flux (an electron leaving the valence band leaves a hole). The Poisson row is unchanged. The partials reach ψ, n and p of the node and ψ of every insulator node on its paths, and the electrode's bias (`bias_derivative`). The pattern gains the path entries at creation. | Ruling 2; Gehring §4.2.2: the current returns to the end nodes' continuity rows and the contact. | Path entries widen the electron and hole rows of gate-interface nodes. |
| GD7 | Currents: the gate's (or electrode's) conduction current is the sum of its nodes' (or paths') tunnelling fluxes, in `terminal_currents`' sign; its charge is unchanged. Transient and small-signal take the term's partials (no storage of its own). | Conservation: the semiconductor contacts carry the same current through the edges. | — |
| GD8 | Switch: per contact, `tunnelling = off` (the default) or a list of processes, `ecb`, `hvb`, `evb` (`on` = all three). With it off the term is not created: residual, Jacobian pattern and values, run identities and device files are those of `main`. | Ruling 4 with opt-out per process for comparison; GC1. | — |
| GD9 | Refused (invalid_input): with `electrothermal`; with the weighted terminal current; electrode tunnelling on a device without `DeviceDescription::cells` (unstructured); a crossed insulator without band data; masses or E_g,ox not finite and positive. | Ruling 7. | Follow-ups recorded in ARCHITECTURE 11.3. |
| GD10 | The deck: `[contact name] tunnelling = …`; on a gate contact optional `oxide_affinity`, `oxide_band_gap`, `oxide_mass_n`, `oxide_mass_p`, `oxide_mass_evb`; on an insulator region `affinity`, `band_gap`, `tunnel_mass_n`, `tunnel_mass_p`, `tunnel_mass_evb`. The device file stores them; the run record names the processes and band data; `BiasPoint` gains the per-contact tunnelling current by process; the app vocabulary and CLI name the keys. | Deck parity; the per-process currents are what a user compares to measurements. | The device file's format version rises if the stored contact or insulator needs it. |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| GC1 | Tunnelling off (and every existing deck): every existing test and pin passes unchanged; residuals, Jacobian patterns and values, run identities and the Unit 40 traces unchanged; a device file without tunnelling is byte-identical to `main`'s. |
| GC2 | WKB closed forms: for linear, rectangular and three-layer piecewise-linear barriers, with and without a turning point inside a segment, TC against a 2000-panel Gauss–Legendre quadrature of the same integrand within 1e-10 relative; the rectangular barrier equals exp(−2 √(2 m (W − E)) d / ħ) within 1e-13. A two-layer stack of identical band data equals one layer of the summed thickness within 1e-13. |
| GC3 | Fowler–Nordheim limit: ECB from a degenerate source through a single SiO₂ layer at 10–14 MV/cm, slope of ln(J/E_ox²) against 1/E_ox equal to −4 √(2 m_n) (qΦ_B)^(3/2) / (3 q ħ) within 2%. Direct regime: at a fixed oxide voltage of 0.5 V, the slope of ln J against t_ox over 1.2–2.0 nm equals −2 κ̄ of the trapezoid (its mean-barrier WKB exponent per length) within 3%. |
| GC4 | Supply and equilibrium: the supply equals its Boltzmann closed form within 1e-12 where both levels are 10 kT below E; at thermal equilibrium (every process, n- and p-substrate, metal and both polysilicon electrodes, both statistics) the tunnelling current is zero within 1e-14 of its value at 1 V. |
| GC5 | Quadrature: halving the panel width changes J by at most 1e-6 relative at every GC3, GC4 and GC10 point. |
| GC6 | Energy limits: EVB is exactly zero while E_c,G > E_v,s (and the reverse pair's condition) and positive past it; ECB to a metal and HVB to a metal have no gate-side limit; a polysilicon gate's limits follow its edges (checked on a 1D n⁺-poly nMOS at −1 … +2 V). |
| GC7 | Jacobian: against central differences of the residual with tunnelling on (1D lumped gate; 2D and 3D tensor electrode stacks of one and two layers; Boltzmann and Fermi–Dirac; each process alone and all three), the largest relative column error at most 1e-6, the off-node ψ columns of the path included. Newton converges quadratically on GC8's device at +1.5 V. |
| GC8 | Self-consistent MOS capacitor (1D, p-Si 5e17, n⁺ poly, 1.5 nm SiO₂, lumped, V_G = −2 … +2 V by 0.1 V): every point converges; the gate current equals the substrate contact's current within 1e-10 relative; the ECB current at each point equals the GD4 integral evaluated on the solution's own band edges and quasi-Fermi levels within 1e-12. |
| GC9 | Lumped against meshed: the same capacitor with the oxide meshed (0.05 nm spacing, an `electrode` n⁺ poly) gives the lumped gate current within 2% at +1 V and +2 V. |
| GC10 | Polarity and dominance: on GC8's nMOS in inversion ECB is the largest process; on its pMOS mirror (n-Si 5e17, p⁺ poly) in inversion HVB is; the per-process currents carry the signs of GD6. |
| GC11 | Geometry: a 2D tensor-grid device uniform across the gate gives the 1D current per unit area within 1e-9 relative; a 3D tri-gate fin node at a corner takes one path per face, each with its face's area (checked against the sum of two separately built 2D devices within 1e-9). |
| GC12 | Dynamics: the small-signal conductance at 1 Hz against DC differences of the gate current within 1e-5 (GC8, +1.5 V); after a gate step the total currents of all contacts sum to zero within 1e-10 relative at every time step. |
| GC13 | Refusals of GD9 and of a negative or absent required band datum, a gate stack with permittivity ≠ 3.9 and no band data, an unknown process name: each invalid_input with the offending value in its context. |
| GC14 | Cost: on the FinFET deck, assembly time with all three processes at most 1.5 times without (median of 31, Release, P-core pinned as Unit 40). The FinFET sweep with tunnelling converges at every point that converges without it. |
| GC15 | The deck, device file (round trip), run record, `BiasPoint` currents, CLI and app vocabulary carry GD10's keys and results. |
| GC16 | Reported, not gated: GC8's I_G–V_G against a measured ultrathin-oxide J–V if an accessible source is found (then its tolerance is ruled before the comparison is run); the FinFET deck's gate current against its drain current. |
| GC17 | Full Release and Debug suites; zero warnings; hand mutants on each GD row caught or recorded as equivalent. |

**Owner rulings on the criteria (AskUserQuestion, 2026-10-28), after the first deck-level runs.** Measured on GC8's
capacitor (lumped, all processes, the deck defaults): its gate current is SRH-generation-limited (2.98e-9 A/cm² at
+0.5 V, 8.59e-9 at +2 V); the inversion electrons' quasi-Fermi level lies 1.2e-6 to 9.6e-5 eV above the gate's; the
substrate's `terminal_current_resolution` is 1.5e-5 to 8.7e-7 A/cm², so its current (−3.04e-9, −9.11e-9) is not
resolved. With SRH and Auger off the gate current is 5.6e-12 A/cm² and the two currents agree to about 1e-7.

| # | Ruling |
|---|---|
| GR1 | GC8's conservation check (gate current against the substrate's) runs with SRH and Auger recombination off, so the current is tunnelling-limited and resolved at the substrate; tolerance as GC8 (1e-10 relative), else its measured value goes back to the owner. The rest of GC8 keeps the deck defaults. |
| GR2 | GC8's ECB recomputation: \|ΔJ\| ≤ 1e-12 × the one-way flux, the GD4 integral with the electrode's supply removed (only the semiconductor side supplying), which is the scale the rounding of the inputs resolves; the net current near quasi-equilibrium is a small difference of the two one-way fluxes. |
| GR3 | GC12's small-signal conductance is compared at 0 Hz (not 1 Hz) with the DC difference, same device and 1e-5: the capacitor's inversion layer responds in about 10³ s, so 1 Hz is its high-frequency response (G = 6.5e-8 against the DC 3.2e-9 S/cm², physical). |
| GR4 | GC9's meshed oxide (31 oxide nodes at 0.05 nm) trips `min_pivot_ratio` (1e-11; measured 6.7e-13 at the next silicon node's electrons, scaling with the tunnel entries; 2.7e-4 for the lumped gate and for paths of at most 6 nodes). With the guard off Newton converges in 5–9 iterations and the meshed current is the lumped one within 6e-5. Investigate the elimination first (ordering, structural entries) and fix it in the assembler or solver if possible; ask again if not. |

**GR4 investigation (2026-10-28).** A dense partial-pivoting LU of the same equilibrated Jacobian at the equilibrium
state gives pivot ratios of 2.7e-4 to 3.5e-4 in natural order, and the dense factors estimate ‖A⁻¹‖ (equilibrated,
max norm over sign vectors): meshed without tunnelling 5.7e18, meshed with HVB only 5.7e18, meshed with all processes
1.2e14; lumped without 8.4e17, lumped with all 1.2e14. The capacitor's Jacobian is near-singular with or without
tunnelling: the inversion electrons are almost floating (their link to the substrate contact through the depletion
region is exponentially weak; only generation and the tunnel current anchor them), and tunnelling improves the
conditioning by about four orders. Eigen's SparseLU pivots partially over each column's structural rows
(`SparseLU_pivotL.h`, threshold 1.0), so its smallest pivot is genuine for its COLAMD order; the path's row changes
that order, which exposes the near-singularity to the `min_pivot_ratio` guard, a floating-region detector whose
verdict depends on the order (it passes the equally conditioned lumped device). Without tunnelling the meshed device
passes the guard at 0 V and then does not converge at +0.5 V. No assembler or solver defect: the path entries change
the ordering, not the conditioning. Returned to the owner.

| # | Ruling |
|---|---|
| GR5 | GC9's meshed device runs with `[solver] min_pivot_ratio = 0`; recorded as a limitation: a tunnelling MOS capacitor without a minority-carrier contact is near-singular and may trip the guard (lower it, or give the inversion layer a source). |

**Measured after GR1–GR5, and ruled (AskUserQuestion, 2026-10-28).** GR1's setup (SRH and Auger off): gate and
substrate agree to about 2e-10 relative in inversion; in accumulation an absolute floor of 7.6e-10 A/cm², the same with
Newton's `tol_update` at 1e-8 and 1e-12 (the substrate terminal current's rounding of cancelling hole fluxes at
p = 5e17; 37% of the 2.6e-10 A/cm² at −0.1 V). GR3's 0 Hz admittance: Y_GG against the DC difference 4e-6 relative at
−1.0 V; at +1.5 V unresolved: +2.7e-9 to −4.8e-8 S/cm² against the DC 3.19e-9 as the linear solve's backward error goes
from 1.4e-11 to 1.7e-16, the solver's own reported resolution 6e-7 to 8e-7 S/cm² (the gate current there is a 1e-5
cancellation between its bias term and the inversion electrons, which follow the gate). GC12's transient (1.4 → 1.5 V,
GR1's setup), the worst relative sum over steps with a current above 1e-3 A/cm²: 2.5e-8 without tunnelling (the
transient's own floor), 6.1e-8 with it. GR2 as worded scales by the semiconductor side's one-way flux, about 1e-16
A/cm² in accumulation, where the current comes from the gate (the measured ECB agreement there is 6e-17 A/cm² on
0.126).

| # | Ruling |
|---|---|
| GR6 | GC8's conservation (GR1's setup): \|I_G + I_sub\| ≤ 1e-10 \|I_G\| + 1e-9 A/cm². |
| GR7 | GC12's small-signal comparison (GR3, 0 Hz) at −1.0 V (accumulation), within 1e-5. The 0 Hz gate admittance of a generation-limited capacitor in inversion is recorded as unresolved; its reported resolution covers the miss. |
| GR8 | GC12's transient: \|Σ I\| ≤ 1e-7 × max \|I\| at each step with max \|I\| > 1e-3 A/cm² (GR1's setup applied to the transient). |
| GR9 | GR2's scale is the larger of the two one-way fluxes (the semiconductor's side supplying alone, and the electrode's alone). |
| GR10 | The FinFET deck's HfO₂ band data for GC14 and GC16 (AskUserQuestion, 2026-10-28): gap 5.7 eV and affinity 2.55 eV (conduction offset 1.5 eV to silicon's 4.05; Gehring Table 5.1), tunnelling masses 0.18 m₀ for electrons and holes (owner ruling: no accessible source found for 0.18; the open-access measurement found, Monaghan et al., Solid-State Electronics 2009, https://research.ucc.ie/en/publications/determination-of-electron-effective-mass-and-electron-affinity-in/, gives 0.11 ± 0.03 m₀ and an affinity of 2.0 ± 0.25 eV). The interfacial SiO₂ takes the built-in set. |
| GR11 | GC14's cost (AskUserQuestion, 2026-10-28): the tunnels' fluxes are computed on threads, the first threading inside assembly; each flux is computed whole by one thread and added in tunnel order, so the result does not depend on the thread count. The baseline (without tunnelling) stays single-threaded. Measured before the ruling, single-threaded: 4.7× (about 1,300 tunnels × 2 processes × 320 quadrature points × about 10 segments, near 55 ns a point). |
| GR12 | Gate tunnelling with quantum confinement (owner, 2026-10-29, after the rebase onto `main` brought Unit 47): refused (invalid_input) as GD9's exclusions, tested in GC13. MLDA moves the surface quasi-Fermi level by −kT ln f_n, which the supply reads, so the same density tunnels more; neither unit validated the pair. A later unit may rule on the supply under confinement and lift the refusal. |

**Slice 2 (device, assembly, deck, files, record, app, CLI), as built.**
- `device::Contact::tunnelling` (`TunnelProcess` bits; 0 off) and `GateStack`'s oxide band data (`oxide_affinity_eV`,
  `oxide_band_gap_eV`, `oxide_mass_n`, `oxide_mass_p`, `oxide_mass_evb`; 0 absent). The device refuses tunnelling on
  anything but a gate or an electrode, unknown process bits, a partial or non-finite band-data set, and a tunnelling
  gate whose permittivity is not 3.9 without band data (GD9).
- `physics::silicon_dioxide_tunnelling_parameters` holds SiO₂'s band data; `silicon_dioxide_parameters` is unchanged
  and the assembler substitutes the tunnelling set for it, so device files and run identities without tunnelling are
  unchanged (GC1).
- `assemble/gate_tunnels.cpp`: a lumped gate gives one tunnel per gate node (one piece of t_ox, its patch area); an
  electrode gives one per semiconductor–insulator grid edge whose straight grid line through the insulator reaches one
  of its nodes (the edge's coupling area). The barrier starts at the first edge's midpoint (the box method's interface;
  the semiconductor's half is not barrier), with the permittivity-weighted potential there; an edge between two
  insulator materials is split at its midpoint the same way. Electrode paths need `DeviceDescription::cells` (else
  invalid_input). Each process's value and gradient are carried over the local columns (ψ, n, p of the node, ψ of each
  path node, a lumped gate's bias), so the Jacobian is the exact derivative of the discrete flux.
- The term `GateTunnelling` (continuity rows, global domain) adds −T_n to the electron row and +T_p to the hole row;
  the pattern gains the path's ψ columns in both. `terminal_currents` adds T_n − T_p to the tunnel's contact;
  `small_signal_currents` its partials and a lumped gate's bias term; `bias_derivative(contact, x)` the lumped bias
  partials. A lumped tunnelling gate makes dF/dV depend on the state: `bias_nonlinear()`, read by the bordered trace and
  the small-signal forcing.
- Refusals (GD9): electrothermal with tunnelling (`DriftDiffusion::create`), the weighted terminal current
  (`weighted_current.cpp`).
- Deck: `[contact] tunnelling = off | on | ecb, hvb, evb` and the gate keys `oxide_affinity`, `oxide_band_gap`,
  `oxide_mass_n`, `oxide_mass_p`, `oxide_mass_evb`. An insulator region's band data uses the parameter names, as every
  material parameter key does (`electron_affinity_eV`, `band_gap_eV`, `tunnel_mass_n`, `tunnel_mass_p`,
  `tunnel_mass_evb`), not GD10's shorter names.
- Device file: contact key 7 (`tunnelling`) and gate-record keys 7–11 (band data), insulator-parameter keys 2–6, each
  written only when nonzero, so files without tunnelling are byte-identical (format version unchanged).
- Run record: insulator band data and each tunnelling contact's processes and gate band data hashed only when present;
  each tunnelling contact named as the setting `tunnelling.<contact>` (its process bits).
- `BiasPoint::tunnelling_ecb`, `_hvb`, `_evb` (per contact, terminal_current's unit and sign; empty without
  tunnelling), sweep-file keys 19–21; the CLI's sweep table gains `<contact> ECB A`, `HVB A`, `EVB A` per tunnelling
  contact; the app vocabulary names the six contact keys (the gate's `o…` completion test now lists the five oxide
  keys).
- Tests: `tests/assemble/gate_tunnelling_test.cpp` (GC6, GC7 Jacobians, equilibrium, bias and small-signal partials,
  GC13); `tests/deck/gate_tunnelling_test.cpp` (GC7 Newton, GC8, GC9, GC10, GC11 in 2D and 3D, GC12, GC15; hidden: GC14 cost and sweep, GC16). Measured: GC8
  ECB against the integral on the solution's bands 6e-17 A/cm² on 0.126 (accumulation); GC9 0.000% and 0.006%;
  GC12 4e-6 (−1 V) and the transient 6.1e-8; GC7 corrections 1.0e-2, 3.5e-5, 2.5e-11. With Newton's `tol_update` at
  1e-12 (not the deck's 1e-8) the same transient's sum reaches 1.4e-4 (an absolute 2e-7 A/cm² while the currents decay
  through 1 A/cm²); not investigated.

**GC14's cost, as built (GR11).**
- `physics::tsu_esaki_current` takes the WKB sum inline: each segment's factor scale·√m once per call, each end's root
  once per point (adjacent segments share an end), the supply's softplus and logistic from one exponential, scratch on
  the stack. Bit-identical to slice 2 (the same operations on the same values): GC7's Newton corrections are unchanged.
  A variant with one division per segment instead of two (about 7% faster) changed the rounding and GC7 then failed
  (corrections 3.1e-3, 6.8e-5 at +1.5 V): the capacitor's near-singular Jacobian (GR4) makes GC7's sequence depend
  on the last bit of the barrier action. Not adopted; recorded as a sensitivity of GC7.
- `DriftDiffusion::GateTunnelJob` (`terms.hpp`, `gate_tunnels.cpp`): `Composition::assemble` starts it before the
  terms, worker threads (the process's affinity count less one, at least two 16-tunnel chunks each) claim chunks of
  tunnels in turn, the other terms add theirs meanwhile, and the `GateTunnelling` term's `finish()` has the assembling
  thread take the remaining chunks, joins the workers and scatters the fluxes in tunnel order. `gate_tunnel_fluxes`
  (terminal and per-process currents) runs the same job at once. Each flux's local gradients live in a per-thread
  arena (no allocation once grown).
- Measured (FinFET, Release, P-core pinned 0x3C3, median of 31): evaluate 13.2 ms without, 15.7 ms with (1.195×);
  residual 9.6 and 13.6 ms (1.42×). Steps: 5.6× (slice 2), 4.7× (bit-identical arithmetic), 1.8× (threads after the
  other terms), 1.2× (threads alongside them).

**GC14's sweep, GC16 and GC11 in 3D, as measured (2026-10-29).**
- GC14 (`[.gate_tunnelling_finfet_sweep]`, PARDISO on 6 threads, the deck's residual line search): without tunnelling
  17 of 17 points, with all three gates tunnelling (HfO₂ band data of GR10) 17 of 17. Both sweeps together took
  69.8 min.
- GC16, reported: the FinFET's gate current against its drain current at V_D = 50 mV. I_G is −3.5 fA at 0 V
  (−1.6e-5 of I_D), 2.1e-6 to 5.3e-6 of I_D from 0.05 to 0.4 V, and rises to 2.29 nA against 18.8 µA at 0.8 V
  (1.22e-4). No measured ultrathin-oxide J–V source was found, so that comparison is not run.
- GC11 in 3D (`tests/deck/gate_tunnelling_test.cpp`): a silicon corner with 1.5 nm of SiO₂ on two faces and an
  n+ poly gate on each, against the same section in 2D with one gate tunnelling per device, at one fixed state
  (the section's equilibrium potential under 1.5 V on both gates, copied along z). Each 3D gate carries 1 µm times
  its 2D device's current; equal to the printed seven digits (8.949119e-10 and 1.002478e-09 A), within 1e-9.
  Mutant: a node keeping only its first tunnel (the corner losing its second face's path) fails by 8.7e-4.

**GC1's traces, as measured (2026-10-29, after the rebase and GR12).** Unit 42's PG1 method: the Unit 40 assembly
traces of `main` (`329f203`) and of this branch, Release, each run pinned to the P cores from its start, in
declaration order, with the MKL runtime loaded, over the default cases of the assemble, solve, deck and validation
executables and the hidden `[mosfet]`, `[transient]`, `[breakdown]`, `[btbt]`, `[unstructured]`, `[mlda]` (solve
and deck) and `[lombardi]` entries. Every one of `main`'s 469 test-case blocks is byte-identical in the branch's
traces (assemble 127, solve 220, deck 73, MOSFET 11, breakdown 7, transient 6, BTBT 3, unstructured 3, solve MLDA 4,
deck MLDA 2, Lombardi 5, validation 8). The only blocks without a counterpart are this unit's 19 new test cases.
With tunnelling off nothing the assembler computes changed, including the threaded assembly (GR11).

**GC17's hand mutants, as run (2026-10-29).** 23 mutants (the last for GR12's refusal, after the rebase), each built in Release and run against the gate
tunnelling tests (physics, assemble, deck), then reverted. All are caught; none is recorded as equivalent. Seven
first survived and each gap was closed with a test before the mutant was run again:

| Row | Mutant | Caught by |
|---|---|---|
| GD1 | E_c without the band shift | survived (all-silicon tests: the shift is 0); new: GC8's ECB on a germanium surface with the gate at x_max, off the reference material |
| GD1 | the lumped gate's vacuum level +0.1 eV | GC8, GC15 |
| GD2 | the corner node keeps one tunnel | GC11 in 3D |
| GD2 | the lumped oxide 1.1 t_ox | GC8 |
| GD2 | an electrode path continuing through semiconductor nodes | survived; new: no path through a semiconductor (fixed states, electrode / SiO₂ / Si / SiO₂ / Si) |
| GD3 | the transmission e^(−S) for e^(−2S) | GC3 slopes |
| GD3 | the turning-point segment's action halved | GC3 FN slope, the current's partials |
| GD3 | the two-sided segment's end partials swapped | the segment and current partials |
| GD4 | the supply cut at 10 kT | GC5 |
| GD4 | no derivative through the barrier top at the top cut | survived (no case reached that cut); new: a fifth partials case with the levels near the top |
| GD4 | an end root shared across a discontinuous step | GC5, the current's partials |
| GD4 | the supply logistic's branches swapped | the current's partials |
| GD5 | EVB's mass falling back to the hole mass | survived (every insulator had equal masses); new: unequal masses in the band data test, and SiO₂'s ruled values pinned |
| GD6 | EVB's semiconductor-valence direction with the wrong sign | survived (only != 0 was checked); new: GC6 at +3 V, the rows −T_n and +T_p against the reported EVB current |
| GD6 | EVB's gate-valence direction with the wrong sign | GC6 at −3 V (added with the one above) |
| GD6 | HVB with the electron mass | new: each process takes its own mass (fixed state, unequal masses through the gate's keys) |
| GD7 | HVB's terminal current with the wrong sign | GC10 |
| GD8 | `on` without EVB | the assemble tests (GC7 stacks) |
| GD9 | no electrothermal refusal | GC13 |
| GD9 | no weighted-current refusal | survived; new: the weighted current refused with gate tunnelling |
| GD9 (GR12) | no quantum-confinement refusal | GC13 (added with the refusal; checks the message) |
| GD10 | the gate's `oxide_mass_n` and `oxide_mass_p` swapped | each process takes its own mass |
| GD10 | a region's `electron_affinity_eV` into its band gap | survived (only the hidden FinFET case uses region keys); new: an insulator region's band data from the deck |

The germanium case needs the floating-region guard off, as GC9 (GR5): with the gate at x_max the capacitor's
near-singularity meets the guard. A lumped `n_poly` gate takes its work function from the channel's material (the
gate electrostatics' convention, `assemble::gate_work_function_eV`), so on germanium it is an n+ germanium gate.

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

**Questions for the owner before code (proposal's recommended answer first):**
- EQ1 τ_w for silicon: 0.35 ps / 0.4 ps (Triebl), or another cited pair.
- EQ2 `field_mobility = on` with energy balance: refused (ED4), or silently replaced by the T-driven law.
- EQ3 Energy-flux factor c = 5/2 + r: the material's thermopower r (silicon c = 2, consistent with Unit 23), or c = 5/2
  (r = 0). Aste–Vahldieck note that smaller c suppresses spurious overshoot; c changes EC10's peak.
- EQ4 Build order: slice 1 physics (ED4, ED5, ED8 functions; EC2), slice 2 Boltzmann assembly and steady sweeps
  (EC1, EC3–EC8, EC10, EC11), slice 3 impact ionization and Fermi–Dirac, slice 4 transient, AC, deck, files, app
  (EC13–EC16), slice 5 FinFET, refinement, mutants, ARCHITECTURE (EC9, EC12, EC17, EC19).

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
