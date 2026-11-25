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
- **Amendment (owner, 2026-11-19: "go with option 1").** The fallback above is now in use.
  - **The breakage.** Runner image `windows-2025-vs2026` 20261116.315 moved Visual Studio from 18.11.12328 to
    18.12.12405. 18.12's default toolset (`VC.Tools.x86.x64`) is 14.52. 14.51 had only ever been 18.11's default,
    never a side-by-side component, so it went with no note in the image's release notes. Since 2026-11-18 about
    17:38 UTC, every job has failed at vcvars: "Toolset directory for version '14.51' was not found".
  - **The fix.** The workflow's "Ensure the MSVC 14.51 toolset" step adds
    `Microsoft.VisualStudio.Component.VC.14.51.x86.x64` with the Visual Studio installer when the image lacks it.
  - **The pin moves from 19.51.36260 to 19.51.36263.** That is the build the 14.51 component installs now. 14.51 was
    serviced in place: the owner's 18.12 has it in the same `14.51.36231` directory, and every local build and test
    since then has used it. The pin stays on 14.51 because 14.52.37025 miscompiles `= {}` (Q11).
  - **Both builds accepted (owner, 2026-11-19).** GitHub rolls images out gradually. PR #94's first release job
    landed on an 18.11 runner, whose default 14.51 is still 19.51.36260, and failed the 36263-only check. The
    assertion now accepts exactly 19.51.36260 or 19.51.36263, nothing else. 36260 is dropped once every runner is on
    18.12.

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
- **Unit 17 — intelligent adaptive mesh refinement and state transfer** ([record](docs/decisions/unit-17.md); `main`
  `72377ad`, PR #81).
  - **AD1–AD10, AC1–AC21, rulings AQ1–AQ5, ID/AR rows:** goal-driven refinement through the generators' size field
    (tensor, 2D and swept meshes), physics indicators, adjoint goal weighting (a transposed solve), one mesh per sweep,
    state transfer, a self-reported certificate; no coarsening or anisotropy.
  - Slice 7 suites: Release 217/217; Debug 198/199 (the known `app.live_check` flake under `-j 6`).
  - **Accepted as limitations (ARCHITECTURE.md 11.3):** AC14, AC15, AC16 (well MOSFET tolerance mode, 1.92× the
    DOFs; tensor budget mode, a tie within AR9's uncertainty, `[!shouldfail]`), and AC8's FinFET target.
- **Unit 17b — cross-run learning for adaptive refinement** ([record](docs/decisions/unit-17b.md); `main` `dae4568`,
  PR #82).
  - **LB1–LB13 (Unit 17's ID8–ID10, AC19, AC21):** opt-in learned first mesh from adapt history (`[output]
    adapt_history`), gradient-boosted trees plus ridge in a versioned `.nitlearn` model, reproducible and
    device-held-out; 22 hand mutants killed.
  - AC19: recall 0.90, precision 0.72; m7 met (a pass saved). **Accepted as a limitation:** m6's DOF miss (1.15×,
    `[!shouldfail]`), the criterion unchanged.
- **Unit 17c — adaptation for derived goals: V_T, g_m and swing** ([record](docs/decisions/unit-17c.md); `main`
  `9fd9fa1`, PR #85).
  - **DG1–DG12, DC1–DC10, DQ1–DQ4, DR rows:** `solve/adapt_derived` with each derived goal's own adjoint over several
    points; the extractor formulas restated in solve and held bit-equal to `analysis/dc.hpp` by a test; V_T absolute
    1 mV, p-channel by −1 scaling.
  - DC10: Release 222/222. DC4 misses for the swing on both halves (its error is in y, DR6). The two derived-goal deck
    cases run in Release only (`deck.adapt_derived_release`; owner, for the Debug CI budget).

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
- **Unit 26d — GUI quick wins** ([record](docs/decisions/unit-26d.md); `main` `9fa3759`, PR #86).
  - **ZD1–ZD8, ZQ1–ZQ4, ZG1–ZG10:** Apply removed (Enter or focus loss), bias/time point labels, colour range
    min/max/auto, drag and drop, taskbar progress and flash, dark title bar, units in the metrics strip, viridis as
    the app's 3D default (`nitcad plot` keeps turbo).
  - App tests 66/66, window tests 24/24, 13 hand mutants caught; the Debug suite left to CI.
- **Unit 57 — generic fields and equations, completed** ([record](docs/decisions/unit-57.md); `main` `2306228`, PR
  #87; accepted by the owner 2026-11-18).
  - **GE1–GE7, GG1–GG9, GQ1–GQ6 (a):** model data in `SystemData`, state in `SystemState`; `DriftDiffusion` retired
    for `assemble::System` and topic headers; derived and edge fields (`derived.hpp`); equilibrium Poisson from the
    same terms (`EquilibriumPoisson` removed); GG8's test-only field χ joins a composition with no change to
    `System` (DofMap and pattern grow by one field, FD, quadratic Newton, k = 0 identity, second-order MMS).
  - Every slice bit-identical (81/81 traced entries at GG8, 38 deck files), FG4 without slowdown, 233/233 Release,
    9/9 GG8 mutants. MSVC 14.52.37025 miscompiles (`= {}` dropped at /O1, /O2): build stays on 14.51 (Q11).
- **Unit 59 — a natural-test norm for reverse bias** ([record](docs/decisions/unit-59.md); `main` `38163aa`, PR
  #101; closed 2026-11-20).
  - **ND1–ND5, NG1–NG7, NQ1–NQ7 (recommended):** opt-in `[solver] natural_norm = weighted` (default
    `convergence`, Unit 30's test unchanged); `System::natural_size` measures both corrections at the iterate,
    held unknowns skipped, densities in asinh with τ = 1e-8 (the only slice 0 candidate to pass NG3 and NG4).
  - Natural search: diode 0.6/0.8/1.0/−5 V in 37/44/50/61 iterations (−20 V still fails, reported); MOSFET
    benchmark 10/10 in 93 (2/10 with Unit 30's norm, a drop bisected to Unit 35's defaults). NG1 86/89 traces
    byte-identical, 3 with added cases only; NG2 38/38 deck files; 8/8 mutants. Suites from CI on `38163aa`
    (Release 248/248, Debug 227/227) and Unit 58's later local runs.
- **Unit 60 — release version policy (PD5) and cross-machine comparison** ([record](docs/decisions/unit-60.md);
  `main` `9093cc8`, PR #96; closed 2026-11-20).
  - **RV/XC/XG/VQ/XQ (recommended):** SemVer 0.1.0 (untagged; only the owner tags), a `v*` tag check in CI;
    `nitcad compare` as a tolerance comparison, not a bit digest (6.8), exit 4 on `differ`; references for
    `pn_diode_1d` and `mosfet_2d` under `tests/cross_machine`.
  - XG4: Release bit-identical across Intel Core Ultra 5 225F and the CI runner's AMD EPYC 7763; Debug against the
    Release references at round-off (largest ratio 2.1e-6, the diode's terminal current). 8/8 mutants.

---

## Open units (records kept here in full)

Merged to `main`, records not yet archived: Units 26e–26i, 58, 61 and 62.

### Unit 26e — hover readout and one default colormap (proposed 2026-11-18; on `app/hover-readout`)

**Owner request (2026-11-18):** "start Hover readout of values; choose viridis or turbo as the default colour map
(GUI)". The hover readout was asked for on 2026-11-16 ("hovering over the 1D/2D/3D model shows no value") and
deferred then. The colormap question left open after 26d is that the app defaults to viridis (ZD8) while `nitcad plot`
keeps turbo, so the same options give two images. The branch starts from `main` `9225e6f`.

**Standing conditions:** as for 26d. The app and `nitcad plot`'s default only, with no change to numerics, results,
files or run identities. Win32 + D3D12 + Direct2D only. Criteria are fixed before code. Release and Debug suites with
`/W4 /WX`, plus hand mutants. WG1 and the self-tests stay as ruled. Nothing is pushed or merged without authorization.

**Today:** a click in Probe mode reports the nearest node's value (`probe_at`), in the inspector and the status bar.
Hovering puts only the position in the status bar's cursor part (`viewport_hover`), on field scenes only. Curve plots
(the Curves tab and the 1D device views) are not pickable: their frame is laid out inside `Renderer::render_plot`.

**Design:**

| # | Change | Where | Behaviour |
|---|---|---|---|
| HD1 | **Hover readout on field views.** Moving the mouse over a pickable field scene (`pickable()`: 2D and 3D field, slice and iso views; not exploded) shows a readout beside the cursor: `<symbol> = <value> <unit>`, then the position in µm (z in 3D) and the node. It is the node Probe would pick: `Prober::pick` with the field's values and hover's 8 px reach. The status bar's cursor part gains the value. | `viewport_hover`; new pure `hover_text` in `probe` | Off the device (no node within reach) the readout hides. Numbers use `number()` as the probe does; a missing value reads "no value". |
| HD2 | **Hover readout on curves.** On a plot scene, the readout shows the data point nearest the cursor in screen distance, over every series, within 12 px: the series name (if the plot has more than one), then `<x label> = <x>` and `<y label> = <y>`. | New `render::plot_frame` and `render::plot_position` (the layout `render_plot` uses, moved out of it, unchanged); new pure `pick_curve` in `probe` | Points that do not lie on the axes (non-positive on a log axis, non-finite) are skipped. The run monitor's live plot is not covered. |
| HD3 | **The readout window.** One Win32 tracking tooltip (`TOOLTIPS_CLASS`, `TTF_TRACK \| TTF_ABSOLUTE`), owned by the main window, placed 16 px right of and below the picked point (kept on screen). It hides on leaving the view (`TrackMouseEvent` / `WM_MOUSELEAVE`, a new `Host::viewport_leave`), while dragging, and when the scene, tab or field changes. **View ▸ Display hover readout** (mnemonic Y, the free letter) turns it on or off. It is on by default and saved in the app state. | `viewport.cpp`, `main_window.cpp`, `app_state` | The click probe, Measure and Cut line are unchanged. The tooltip never takes focus or clicks. |
| HD4 | **One default colormap: viridis.** `nitcad plot`'s default `--colormap` becomes viridis, the app's since 26d, so "same options, same image" holds with no options. turbo, cividis, coolwarm and soft stay selectable. The difference view keeps cool-warm unless one is chosen (its sentinel becomes the default, viridis; turbo is then a choice there). | `cli/plot.cpp` option default and XD6 sentinel; `view_settings.cpp` (`takes_coolwarm`, `plot_arguments` omit the default); `self_test.cpp` | The render layer's own defaults (`Colorbar::map`) and its reference images are unchanged: they name their maps. A `nitcad plot` run without `--colormap` draws a different PNG from before. This is the user-visible change. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| HQ1 | Where the readout shows | (a) A tooltip at the cursor, and the status bar. (b) The status bar only. (c) A label drawn into the scene. | **(a)**: at the eye's position. (c) would put an app overlay into the renderer and `nitcad plot`'s drawing path. |
| HQ2 | Which value on a field view | (a) The nearest node's value, as Probe gives. (b) Interpolated in the cell under the cursor. | **(a)**: the same number a click gives. 3D has no cell locator (M16), and nodes are the solved values. |
| HQ3 | Which point on a curve | (a) The nearest data point of any series. (b) Every series interpolated at the cursor's x. | **(a)**: solved points only, no interpolated values that were not computed. |
| HQ4 | A switch for it | (a) A View menu switch, on by default, remembered. (b) Always on. | **(a)**: a tooltip that follows the mouse gets in the way of screenshots. |
| HQ5 | The default colormap | (a) Viridis in the app and `nitcad plot`. (b) Keep the split (app viridis, `nitcad plot` turbo). (c) Turbo in both (undo ZD8). | **(a)**: the owner chose viridis on 2026-11-15. It is perceptually uniform and colour-blind safe. It ends the split that makes 26d's difference-view and equivalent-command special cases. |
| HQ6 | A marker on the hovered curve point | (a) None: the readout names the point and sits by it. (b) A ring drawn on the point. | **(a)**: (b) needs a scene overlay in the renderer. It could be a later item. |

**Owner ruling (2026-11-18):** ruled under the standing instruction "always choose recommended" (2026-11-11): HQ1 (a),
HQ2 (a), HQ3 (a), HQ4 (a), HQ5 (a), HQ6 (a). Implementation authorized by the request ("start").

**Acceptance criteria (frozen before code; tests in `tests/app` and `tests/render`, hidden windows as today):**

| # | Criterion |
|---|---|
| HG1 | `hover_text`: a 2D hit (symbol, value, unit, x, y, node; no z), a 3D hit (with z), a NaN value ("no value"), and a field without a unit, each against expected strings. |
| HG2 | `plot_frame`/`plot_position`: for a known scene and size, the frame equals today's layout (with and without a title), and the corners map to the axis limits on linear and log axes. `pick_curve`: nearest across two series; a log axis with a non-positive point skipped; non-finite skipped; nothing beyond 12 px; one series gives no series line. |
| HG3 | The move of the plot layout out of `render_plot` is bit-identical: every render reference image (`band_diagram`, `iv_plot`, the cut-band plots) matches unchanged. |
| HG4 | Window tests: `WM_MOUSEMOVE` over a 2D result's device shows the tooltip with HG1's text for the node Probe picks there; over empty space it hides; `WM_MOUSELEAVE` hides it; with the switch off nothing shows; during a drag nothing shows; on the Curves tab a move onto a curve point shows its HD2 text; a tab change hides it. The setting survives a save and load of the app state. |
| HG5 | Colormap: `nitcad plot` with no `--colormap` gives the same PNG as with `--colormap viridis`; `--colormap turbo` still gives turbo; the difference view with no colormap is cool-warm and with turbo chosen is turbo. A default app field view's equivalent command names no colormap, and run through `nitcad plot` it gives the app's image (WG4 self-test default). |
| HG6 | No change elsewhere: every render reference passes unchanged; app tests change only where HD1–HD4 change what they check; numerics untouched (no solver, physics or io source changed); WG1 self-tests within their ruled bounds. |
| HG7 | Hover costs no more than a frame: on the largest bundled 3D example's surface view, the median `viewport_hover` (pick + text + tooltip update) over 200 positions is under 16 ms in Release. The time is recorded. |
| HG8 | Hand mutants, each killed by a test: a hover pick ignoring the values (old behaviour), z shown in 2D, the curve pick not skipping log non-positives, a 20 px reach, the layout's margin changed in `plot_frame`, the tooltip not hidden on leave, the switch ignored, the switch not saved, the difference sentinel kept at turbo, the `nitcad plot` default kept at turbo. |

**Slices (each one local commit):** 26e·1 HD2's render layout and the pure `hover_text` / `pick_curve` (HG1–HG3);
26e·2 HD1 + HD3 in the window (HG4, HG7); 26e·3 HD4 (HG5); 26e·4 suites, mutants (HG6, HG8), the ARCHITECTURE "as
built" entry.

**Implementation (2026-11-18).** Slices: 26e·1 `bdd3af9`, 26e·2 `b72a550`, 26e·3 `fad0b0a`, 26e·4 this record and
ARCHITECTURE's "as built" entry (6.14).

Choices made in building, within the design:
- **HD3:** the menu item is **View ▸ Display hover readout** (mnemonic Y). Every letter of "Hover readout" was already
  a View mnemonic, and the menu test requires them to be unique. A hidden window (the tests) keeps the readout's
  state and text but never activates the tooltip, so test runs put nothing on the desktop. A view with no field (the
  structure, the mesh) shows the position and node only.
- **HD4:** the difference view's sentinel is the default, viridis, in both the app (`takes_coolwarm`) and `nitcad
  plot` (XD6). Choosing viridis there therefore still shows cool-warm, as turbo did before 26e.
- **Tests changed (HG6):** the `plot_arguments` expectations drop `--colormap viridis` (now the default; one case
  names cividis to keep the option covered); 26d's ZG10 case is rewritten for one default; the self-test's third
  view names turbo instead of viridis. The 26d film tests that name `--colormap viridis` stay as they are.
- **Found and fixed (owner, 2026-11-18: "fix it"):** a 1D diode sweep from 0 V had no I-V plot with the cathode
  sensed ("a log axis needs a positive value"). Its current is exactly 0 at 0 V and negative after that, and
  `render::curve_plot` drew a log-axis curve by its magnitude only when every value was strictly negative. Now it
  does so when the curve has no positive value and at least one negative one; the zeros stay off the axis, as
  non-positive values always were. Mixed-sign curves and linear axes are unchanged, and no curve that drew before
  draws differently (the case was an error). Tests: `render scene: plots fit their axes…` (from 0 V, mixed, linear)
  and `app views: a diode's I-V from 0 V on a log axis` (both contacts); the old condition fails both. The hover
  window test's sweep starts at 0 V again.

Results: Release app tests 71/71 test cases, window tests 27/27 (24 + 3 new), render 78/78 with every reference image
unchanged (HG3), `nitcad plot` deck tests 12/12. HG7: median hover 0.067 ms, maximum 0.78 ms over 200 positions on
`examples/mosfet_3d.deck`'s doping surface (1,911 nodes). HG8: 10 hand mutants, all killed (the curve pick's own
non-finite skip is redundant, since a NaN position fails the reach test, so that mutant is placed in
`axis_fraction` instead). A first mutant round was discarded: its restores kept old file times, so the next builds
kept earlier mutants. It was rerun with each file touched after restore. **Full Release ctest not run** (owner,
2026-11-18: "full suite is not required rn"), and the Debug suite is left to CI. HG6 therefore stands on the suites
above, not on the whole tree.

### Unit 26f — enhanced 3D look (proposed 2026-11-18; on `app/visual-prototype`)

**Owner request (2026-11-18):** "improve the visuals of 3D model ... the GUI of 3D model looks too basic". Of the
eight suggestions, the owner asked for a prototype of the first three: lighting, feature edges and a gradient
background. Prototype `ff02613` sits behind the environment variable `NITCAD_VISUAL_PROTO=1`. On its FinFET and
nmosfet screenshots the owner ruled "Its good" and asked for this proposal. The branch is stacked on
`app/hover-readout` (26e, not yet merged), so it rebases onto `main` once 26e merges.

**Standing conditions:** as for 26e. The drawing of 3D views only, with no change to numerics, results, files or run
identities. Win32 + D3D12 + Direct2D only. Criteria are fixed before code. Release and Debug suites with `/W4 /WX`,
plus hand mutants. WG1 and the self-tests stay as ruled. Nothing is pushed or merged without authorization.

**Today (`main` and 26e):** 3D triangles get one directional light from over the viewer's shoulder, lit two-sided as
`ambient 0.38 + 0.62·|n·l|`, on a white clear. Region solids show no edges unless **Mesh edges** is on. The prototype
adds the following on 3D scenes only:
- key, fill and rim lights, a sky/ground hemisphere ambient and a soft highlight (reduced on field maps);
- dark depth-tested feature edges where faces turn by more than 30°, and at each batch's open boundary;
- a vertical gradient over the view rectangle.

It has three shortcuts: the edges are recomputed on every frame, the eye direction is `-forward` even in
perspective, and the switch is an environment variable.

**Design:**

| # | Change | Where | Behaviour |
|---|---|---|---|
| KD1 | **A look on the scene.** `enum class Look { classic, enhanced }` and `FieldScene::look`. The render layer's default stays **classic**. The renderer branches on the scene's look, and the environment variable is removed. | `render/scene.hpp`, `renderer.cpp` | Classic draws bit-identically to today. 2D scenes ignore the look. |
| KD2 | **Enhanced lighting.** The prototype's model: key (today's light), fill (lower right), rim `pow(1 − n·v, 3)`, hemisphere ambient from the camera's up, and a Blinn–Phong highlight (0.22 on solids, 0.08 on field maps). Two-sided by turning the normal toward the eye. **Per-pixel eye direction:** the mesh shader passes world position, and the constants carry the camera position (perspective) or `-forward` (orthographic). | `scene.hlsl`; a new pure `frame_lighting(scene, view)` in the renderer, which fills those constants | The colormap is sampled exactly as today, and the lighting only scales it, as classic does. The hover readout and Probe still report the solved values. |
| KD3 | **Feature edges, computed once.** A new pure `render::feature_edges(const TriangleBatch&, float angle_deg)` returns segments: edges shared by two faces whose normals differ by more than the angle (either winding), and edges of only one face (or of more than two). The output is sorted, so images are deterministic. `add_feature_edges(FieldScene&)` appends one `LineBatch "feature edges"` (ink 0.08/0.09/0.11, alpha 0.85, 1.7 px, depth-tested) built from the region solids, field surfaces, slices and electrodes. Arrows, streamlines, isosurfaces and the outline shell are skipped. | `render/geometry` or a new `render/edges.hpp/.cpp`; the 3D scene builders call it when `look == enhanced` | Built with the scene, not per frame: a camera move redraws without recomputing. Clipping and hidden regions rebuild the scene, and with it the edges. |
| KD4 | **Gradient background.** In enhanced 3D, the view rectangle is filled before the scene with a vertical gradient (0.86/0.89/0.93 at the top to 0.99/0.99/1.0 at the bottom). The title band and colorbar panel stay white. | `background_ps` in `scene.hlsl`, one PSO | 2D views and plots are unchanged. |
| KD5 | **The app setting.** An inspector check box **Enhanced 3D look**, on by default, saved in the app state with the other view settings. Off gives classic. The equivalent command adds `--look classic` only when it is off. | `view_settings`, `app_state`, `main_window.cpp`, `self_test.cpp` | Shown on 3D views only (as **Perspective** is). |
| KD6 | **`nitcad plot --look classic\|enhanced`**, default **enhanced**, so the app and `nitcad plot` give the same image from the same options (as HD4 did for the colormap). | `cli/plot.cpp` | A 3D `nitcad plot` run without `--look` draws a different PNG from before. This is the user-visible change. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| KQ1 | The default look | (a) Enhanced in the app and `nitcad plot`. The render layer's default stays classic, existing references unchanged, two new enhanced references. (b) Classic everywhere, enhanced opt-in. (c) Enhanced everywhere, including the render layer, with every 3D reference rebuilt. | **(a)**: users get the look the owner approved. As with HD4, the references name what they draw, so no owner-approved rebuild is needed. (c) churns references with no gain in what they test. |
| KQ2 | The edge rule | (a) Dihedral angle above 30°, plus open batch boundaries. (b) (a) plus view-dependent silhouettes. | **(a)**: device geometry is mostly boxes and extrusions, so (a) draws the outline at any angle. It can also be cached, while silhouettes need work on every frame. |
| KQ3 | Which batches get edges | (a) Region solids, field surfaces, slices, electrodes. (b) Every lit batch. | **(a)**: arrow cones, streamline tubes and marching-tetrahedra isosurfaces have many facets over 30°, and outlining them is clutter. |
| KQ4 | Where the edges are computed | (a) Once, when the scene is built (render layer, pure, testable). (b) Per frame in the renderer (the prototype). | **(a)**: rotating costs nothing extra, and the function can be unit tested. |
| KQ5 | The background | (a) A gradient over the 3D view only. (b) Over the whole image. (c) A dark viewport option as well. | **(a)**: what the owner approved. A dark viewport is a later item, with the dark title bar (26d). |
| KQ6 | Highlights on field maps | (a) Reduced (0.08), as prototyped. (b) None. | **(a)**: as approved. Lighting already scales the colormap in classic. Exact values come from the colorbar, the hover readout and Probe. |

**Owner ruling (2026-11-18):** ruled under the standing instruction "always choose recommended" (2026-11-11): KQ1 (a),
KQ2 (a), KQ3 (a), KQ4 (a), KQ5 (a), KQ6 (a). The request was for the proposal, so implementation waits for the
owner's go-ahead.

**Acceptance criteria (frozen before code; tests in `tests/render` and `tests/app`, hidden windows as today):**

| # | Criterion |
|---|---|
| KG1 | Classic is bit-identical: every existing render reference image matches unchanged, the scene's default look is classic, and 2D scenes draw identically under both looks. |
| KG2 | `feature_edges`: a closed cube of 12 triangles gives exactly its 12 edges. A flat square split into 8 triangles gives only its 4 sides (8 segments). A fold at 20° gives no edge at 30°, and one at 40° gives the fold. A flat strip with mixed winding gives only its boundary. Degenerate triangles are ignored. Output is identical across two calls and with the triangle order shuffled. |
| KG3 | `frame_lighting`: orthographic gives eye `-forward` and a perspective flag of 0; perspective gives the camera position and a flag of 1. Key, fill and up are unit vectors, the fill on the opposite side of the eye from the key. Classic leaves the KD2 constants zero. |
| KG4 | Two new references: `examples/mosfet_3d.deck`'s solid view and potential field view, enhanced, compared under the references' existing rules. Each passes, and each differs from its classic render. |
| KG5 | App and CLI: the check box defaults on and survives a save and load of the app state. Off, the view is classic. The equivalent command omits `--look` when on and gives `--look classic` when off, and run through `nitcad plot` it gives the app's image (the WG4 self-test default). `nitcad plot` with no `--look` equals `--look enhanced`. The check box is hidden on 2D views. |
| KG6 | Cost: building the edges for the largest bundled 3D example's field scene takes under 50 ms in Release, and a camera-move redraw is within 10% of classic's (median of 50 frames). Both times are recorded. |
| KG7 | No change elsewhere: no solver, physics or io source changed; app tests change only where KD5 changes what they check; WG1 self-tests within their ruled bounds. |
| KG8 | Hand mutants, each killed by a test: angle 30° → 60°; the `abs` dropped from the dihedral test (mixed winding); open boundaries not drawn; arrows given edges; output not sorted; the renderer ignoring `look`; the background drawn in 2D; the perspective eye replaced by `-forward`; the setting not saved; the `nitcad plot` default left at classic. |

**Slices (each one local commit):** 26f·1 the render layer: `Look`, `feature_edges` / `add_feature_edges`,
`frame_lighting`, per-pixel eye, background, the environment variable removed (KG1–KG4); 26f·2 the app check box and
`nitcad plot --look` (KG5, KG6); 26f·3 suites, mutants (KG7, KG8), the ARCHITECTURE "as built" entry.

**Implementation authorized (2026-11-18):** "choose recommended and start".

**Results (2026-11-18; slices `fdb7886`, `871c714` and the slice 3 commit; Release, `/W4 /WX`, MSVC 14.51):**

| # | Result |
|---|---|
| KG1 | **Met.** Every existing render reference matches unchanged. The scene's default look is classic. The 2D potential map draws byte-identically under both looks. Outside the suite, four `nitcad plot` classic renders (FinFET and nmosfet decks, solid and field) on the RTX 5060 Ti are byte-identical to the images drawn before the unit. |
| KG2 | **Met.** `tests/render/edges_test.cpp`: the cube's 12 edges; the flat 2 × 2 square's 8 boundary segments; the 20° fold gives no edge and the 40° fold gives one (none at 60°); mixed winding gives only the boundary; degenerate triangles are ignored; output is sorted and identical after shuffling and rotating triangles; −0 and +0 are one position. |
| KG3 | **Met.** `frame_lighting`: classic and 2D give zeros; orthographic gives eye = the view direction within 1e-6, with unit key, fill, up and eye, key from the left, fill from the right; perspective gives the camera's position, outside the bounds along the view direction within 1e-5. |
| KG4 | **Met, with one wording deviation.** The two references are drawn from the image test's own solved 3D MOSFET (the device all its 3D references use), not from `examples/mosfet_3d.deck`. Each passes the references' comparison rules and differs from its classic render in more than 5% of pixels. |
| KG5 | **Met.** The app-views, app-state and window tests: the check box defaults on, is hidden on 2D devices and the Curves tab, and is saved (`enhanced_look`); a new device's views take the saved choice. The equivalent command omits `--look` when on and gives `--look classic` when off; in both cases `nitcad plot` draws the app's image exactly (WARP). `nitcad plot` without `--look` equals `--look enhanced` and differs from classic; `--look shiny` is refused. The self-test accepts `--look`. |
| KG6 | **Met (redraw clause on the hardware adapter, by the owner's ruling (b) below).** On `examples/mosfet_3d.deck`'s doping surface (1,911 nodes, 776 segments), building the edges takes a median of 1.75–1.91 ms (bound 50 ms). A camera-move redraw (median of 50 drag steps in the hidden test window, which draws on WARP): enhanced 13.7–14.2 ms against classic 10.9–11.1 ms, 1.25–1.28× (bound 1.10×). The prototype shader was 1.47×; the cost was cut by computing powers by multiplication and the orthographic half vector once (no pixel change on WARP). What is left is the per-pixel lighting, which WARP runs on the CPU. The edges cost about 0.3 ms and the background about 0.4 ms. The same draw on the RTX 5060 Ti, off screen: 1.02–1.07× classic. After the ruling, the test draws the scene off screen on the hardware adapter with its camera turned 1° before each draw (no rebuild, no vsync): 3.80 ms enhanced against 3.84 ms classic, 0.99× (0.85–1.01× over four runs). It is skipped where there is no hardware adapter (CI), and it still reports the WARP drag (1.22–1.29×). |
| KG7 | **Met** on the render, app and CLI entries (51 of 51, Release); the full Release ctest was stopped at the owner's word ("dont run full suite just push", 2026-11-18) with 229 of 235 entries passed and none failed (`verify_full.tvb04` and five queued entries not finished). The Debug suite is left to CI. No solver, physics or io source changed. App tests changed only where KD5 adds a check box, a setting and a flag. |
| KG8 | **Met.** 10 hand mutants, each killed: angle 30° → 60°; the `abs` dropped (mixed winding); open boundaries not drawn; arrows given edges; output not sorted; the renderer ignoring the look; the enhanced look drawn in 2D; the perspective eye replaced by `-forward`; the setting not saved; the `nitcad plot` default left at classic. Each file was touched after its restore, so no mutant stayed in a later build. |

**Open for the owner (KG6):** the frozen bound is 1.10× and does not name an adapter. Measured on the test
window's WARP adapter it is unmet (1.25–1.28×); on the hardware adapter, which the app uses by default, it is met
(1.02–1.07×). Options: (a) accept KG6's redraw clause as a limitation on WARP, keeping the `[!shouldfail]` test; (b)
rule that the clause is measured on the hardware adapter, and turn the test's check to it (skipped where there is
no GPU, as in CI); (c) cheaper enhanced shading on WARP (for example no highlight or rim there), which changes the
look under WARP. Recommendation: **(a)**. It keeps the criterion as frozen, users draw on the GPU, and (b) would
change which measurement decides the criterion after the results were seen.

**Owner ruling (2026-11-18):** "b". The redraw clause is measured on the hardware adapter, the app's default; the
test checks it there and skips without one. The WARP figure is kept in the record, not as a gate.

### Unit 26g — materials, depth and finish (proposed 2026-11-18; on `app/materials-depth`)

**Owner request (2026-11-18):** "make plans for 4 - 8" of the eight 3D-view suggestions. This unit covers 4
(metallic contacts, glass-like oxide), 5 (ambient occlusion, a soft shadow) and 8 (softer mesh edges, contour
lines, the colourbar's frame). 6 and 7 are Unit 26h. The branch is stacked on `app/visual-prototype` (26f).

**Standing conditions:** as for 26f. Everything here belongs to the **enhanced** look only. Classic, the render
layer's default, draws bit-identically, so every existing reference is unchanged. No change to numerics, results,
files or run identities. Win32 + D3D12 + Direct2D only. Criteria are fixed before code. Release and Debug suites
with `/W4 /WX`, plus hand mutants. Speed is measured on the hardware adapter (KG6 ruling (b)), and WARP figures are
recorded only. Nothing is pushed or merged without authorization.

**Today (26f):**
- Contacts and electrodes are one flat batch, `"electrode"`, coloured 0.33/0.35/0.38.
- The oxide is a translucent flat batch at `insulator_alpha`, drawn by weighted blended OIT.
- Enhanced shading treats both like any solid.
- There is no occlusion or shadow.
- Mesh edges are 0.15/0.17/0.2 at alpha 0.75 (3D), 1 px.
- Contours (`--contours`, off by default) are a dark line batch.
- The colourbar already has a 1 px frame (0.3/0.32/0.35, Unit 25a), so item 8's frame exists. Only its look in
  enhanced is in question (MQ6).

**Design:**

| # | Change | Where | Behaviour |
|---|---|---|---|
| MD1 | **A finish per batch.** `enum class Finish { matte, metal, glass }` and `TriangleBatch::finish` (matte by default). The scene builders mark `"electrode"` metal and the translucent insulators glass. The renderer passes it in the frame flags (bits 3–4). | `render/scene.hpp`, `scene.cpp`, `renderer.cpp` | Classic ignores it. |
| MD2 | **Metal (item 4).** In enhanced, a metal batch's base colour becomes aluminium 0.70/0.72/0.76. Its highlight is stronger and tighter (0.55, `(n·h)⁶⁴`), and it gets a fake environment: the sky/ground ambient taken along the reflected eye ray instead of the normal, so faces turned up read bright and faces turned down read dark, as brushed metal does. | `scene.hlsl` | Contact pads and gate electrodes alike (MQ1). |
| MD3 | **Glass (item 4).** In enhanced, a glass batch's alpha rises toward grazing angles, `a' = a + (1 − a)·0.55·(1 − n·v)³` (Schlick-like), and its colour lifts by the same term, so the oxide's outline shows where it is seen edge-on while face-on it stays as clear as today. | `scene.hlsl` (transparent pass) | Opaque insulators (`--opaque`) stay matte. |
| MD4 | **Ambient occlusion (item 5).** A screen-space pass after the opaque triangles, before the lines and the OIT composite. It reads the depth buffer (D32 made typeless with an SRV), takes 12 samples in a hemisphere at a radius of 2% of the view's height, darkens by up to 35% on solids and 15% on field maps (a mask the opaque pass writes in the colour target's alpha), and is blurred 4 × 4. | `post.hlsl` (`ao_ps`, `ao_blur_ps`), `renderer.cpp`: one R8 target at the supersampled size | Corners where a contact meets the silicon, the steps of the oxide and the gate's foot read as recessed. Translucent glass is not occluded and does not occlude. |
| MD5 | **A ground shadow (item 5).** A soft shadow under the device: a quad 4% of the height below `lo.y`, on the camera's floor plane, its alpha falling off from 0.28 at the footprint's centre to 0 one twelfth of the extent outside it. Drawn before the scene, over the gradient. | `scene.hlsl` (`shadow_ps`), the footprint from the scene's bounds | It follows the camera. It is not drawn when the view looks up from below (elevation < −5°), or on clipped or exploded views. |
| MD6 | **Softer mesh edges and contours (item 8).** In enhanced: 3D mesh edges 0.25/0.28/0.32 at alpha 0.35, 1 px; contours 1.2 px at alpha 0.55, in the ink of the feature edges. | `scene.cpp` (the batch colours and widths read the look through `set_look`) | Contours stay off by default (MQ5). |
| MD7 | **The colourbar in enhanced (item 8).** Its frame is kept and drawn 1 px lighter (0.45/0.47/0.5), with tick marks inside the bar's right edge as well as outside. The panel stays white. | `renderer.cpp` `draw_colorbar` | Classic's colourbar is unchanged. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| MQ1 | The metal's colour | (a) Aluminium grey for every contact and electrode. (b) Gold for contacts and aluminium for gates. (c) By the deck's contact kind and work function. | **(a)**: one finish that reads as metal. Colour codes for contact kinds would be a new legend for no physics gain. |
| MQ2 | Occlusion technique | (a) Screen-space AO (MD4). (b) Precomputed per-vertex AO from the mesh. (c) None, shadow only. | **(a)**: it works for every view, clip and region visibility with no per-scene cost. (b) would need rays through the mesh on every rebuild. |
| MQ3 | AO on field maps | (a) Reduced (15%) through the mask. (b) Full. (c) None on field maps. | **(a)**: the corners still read in depth, and colours shift less than the lighting already shifts them. Values come from the colourbar, hover and Probe. |
| MQ4 | The ground shadow | (a) A soft footprint shadow (MD5). (b) A real shadow map from the key light. (c) None. | **(a)**: it grounds the device at the cost of one quad. A shadow map needs a second depth pass, and its hard shapes would compete with the field. |
| MQ5 | Contours by default (item 8's suggestion) | (a) Keep them off and restyle them (MD6). (b) On by default on enhanced 3D field views, with a new `--no-contours`. | **(a)**: a default overlay changes every field image and needs an off switch in the app and CLI. On log fields with steep doping, contour lines pile up at junctions. It stays one click away. |
| MQ6 | The colourbar | (a) MD7's lighter frame and inner ticks. (b) Unchanged. | **(a)**: it matches the lighter enhanced view at almost no cost. |
| MQ7 | 26f's enhanced references, which 26g changes | (a) Rebuilt in 26g·4, after the owner approves the new images. (b) Kept, with 26g's finishes behind a second switch so the 26f look stays reproducible. | **(a)**: one enhanced look, not two. The 26f look has no users to keep it for, and a second switch would double the test matrix. |

**Owner ruling (2026-11-18):** ruled under the standing instruction "always choose recommended" (2026-11-11): MQ1 (a),
MQ2 (a), MQ3 (a), MQ4 (a), MQ5 (a), MQ6 (a), MQ7 (a) (the rebuild itself still shown to the owner first). The request was for plans, so implementation waits for the owner's
go-ahead.

**Acceptance criteria (frozen before code; tests in `tests/render` and `tests/app`):**

| # | Criterion |
|---|---|
| MG1 | Classic is bit-identical: every classic render reference matches unchanged. Scenes default to the matte finish, and a classic scene with metal and glass batches draws exactly as with matte ones, to the bit. 26f's two enhanced references change by design (MQ7): they are rebuilt in 26g·4 only after the owner has seen the new images and approved them. |
| MG2 | Finishes: `solid_view` and `field_view` mark electrodes metal and translucent insulators glass; opaque insulators, region solids, field batches, slices, arrows and isosurfaces stay matte. |
| MG3 | Metal and glass in images (WARP): on the enhanced 3D solid, the mean colour of the gate electrode's top face is lighter than classic's and lighter than its front face. Along a line across the oxide, the alpha-composited oxide is at least 1.5× as opaque within 3 px of its silhouette as at its centre. |
| MG4 | AO in images (WARP): on the enhanced 3D solid, pixels on the silicon within 4 px of the gate electrode's foot are at least 8% darker with AO than without. A flat face away from any step changes by under 2%. Field-map pixels change by at most half as much as solid pixels at the same occlusion (the mask). Glass pixels are unaffected by AO under them. |
| MG5 | The shadow: present below the device on the default camera (a darker band under `lo.y` than the gradient alone), absent with elevation below −5°, and on clipped and exploded views. |
| MG6 | Three new WARP references (enhanced): the 3D solid (metal, glass, AO, shadow), the 3D potential surface, and the 3D field with mesh edges and contours. Each passes the existing rules. |
| MG7 | Cost on the hardware adapter (off screen, the KG6 method, `examples/mosfet_3d.deck`): the 26g enhanced draw within 1.25× of 26f's enhanced draw. The AO passes alone under 2 ms at 1280 × 960. WARP times recorded only. |
| MG8 | No change elsewhere: no solver, physics or io source changed; app tests change only where the images do; the app and `nitcad plot` still give the same image from the same options (`same_as_plot`, WG4 self-test). |
| MG9 | Hand mutants, each killed by a test: electrodes left matte; the glass term's power 3 → 1; glass applied to opaque insulators; the AO mask ignored (full AO on field maps); AO sampling glass; the AO radius 2% → 10%; the shadow drawn when looking up; metal shading in classic; the mesh-edge alpha left at 0.75 in enhanced; the colourbar change applied in classic. |

**Owner ruling MQ8 (2026-11-18, during 26g·1): the glass on flat faces.** Measured, MD3's grazing term is constant
over each flat face of the oxide box, so on the default view the top face only turns a little paler. There is no rim at
the silhouette, and MG3's glass half (1.5× as opaque within 3 px of the silhouette) cannot be met by flat faces. Asked,
the owner ruled the recommended option, glass rim edges:
- **MD3, amended:** the face term stays as designed. In addition, in enhanced, `set_look` puts the feature edges of
  glass batches into a second line batch, `"glass edges"`: a light rim (0.97/0.98/1.0, alpha 0.9, 1.5 px, depth-tested)
  in place of the dark ink.
- **MG3, glass half, restated:** on the enhanced 3D solid, pixels within 3 px of the oxide's top edges (away from the
  gate) are lighter than the same pixels with the oxide matte. The test asks for a mean luminance more than 10 higher.
- **MG9 gains a mutant:** the glass rim drawn in the dark ink.

**Owner ruling MQ9 (2026-11-18, during 26g·2): AO dropped.** MD4 was built as designed and measured on the test
MOSFET's enhanced 3D solid against the same renderer without it. Its pipeline worked: a forced constant darkening
showed in the image. But the scene gave it almost nothing to occlude:
- The gate's foot sits on the translucent oxide, which writes no depth (by MD4), so in depth the gate is 0.08 µm
  above the silicon. Change: 0%.
- The contact pads are slabs 0.06 µm thick (3% of the extent). Their undersides face away from the camera, and their
  visible faces lie beyond the radius (2% of the view's height, about 0.036 µm here). Change: 0% at 2%, 0.7% at 5%.
- Flat faces: 0%.

Asked, the owner ruled the recommended option: **drop AO from 26g.**
- MD4 and its code are removed. MG4 and the AO half of MG7 are withdrawn, and the two AO mutants (the mask ignored,
  AO sampling glass) and the radius mutant are dropped from MG9.
- Depth comes from the ground shadow (MD5) and 26f's lighting and edges.
- AO may return in a later unit if scenes with real recesses appear (trenches, a FinFET gate wrapping its fin).

**Slices (each one local commit):** 26g·1 finishes, metal and glass (MD1–MD3; MG2, MG3); 26g·2 the shadow, mesh edges,
contours, colourbar (MD5–MD7; MG5, MG6, MG7's whole-draw cost); 26g·3 suites, mutants (MG1, MG8, MG9), the
ARCHITECTURE entry. (Originally four slices, with AO as 26g·2; renumbered after MQ9.)

**Results so far (2026-11-19):**
- 26g·1 `d0b9e55`, MQ9 `62e6378`, 26g·2 `cbd89df`. The owner approved the three enhanced images on 2026-11-18,
  and they were then written as references (`f8faf26`, MQ7, MG6): `mosfet3d_solid_enhanced` and
  `mosfet3d_surface_potential_enhanced` (rebuilt), and `mosfet3d_field_mesh_contours_enhanced` (new). Each is
  byte-identical to the image approved.
- MG1–MG3, MG5 and MG6 are met by the new render tests. MG1's reference half is met by the classic references
  (every one unchanged in the render suite run at 26g·1).
- MG7, on the RTX 5060 Ti by KG6's method: 26g enhanced 5.48 ms against 26f's 5.56 ms (medians of 50), so **0.99×**
  (bound 1.25×). 26f's draw is the same scene with matte batches and no shadow: exactly 26f's triangles and lines. On
  WARP it is 1.08×, recorded only. KG6 still holds at 1.04× classic (bound 1.10×).
- CI fix (PR #92, 2026-11-19). The CI runner's "Microsoft Basic Render Driver" (vendor 0x1414, device 0x8c) is WARP
  behind a display adapter, but it is not flagged as software. `Adapter::hardware` picked it, so the hardware-adapter
  tests ran on WARP and failed: KG6 at 1.59× in Debug and 1.26× in Release. The hardware search now skips it, so those
  tests skip on CI as the KG6 (b) ruling intended. Owner ruling: KG6's edges bound (< 50 ms) is checked in Release only.
  The unoptimised Debug build took 113 ms on CI (Release: 4.9 ms on CI, 1.9 ms on the RTX machine), and Debug now only
  records the time.
- CI fix 2 (PRs #92, #93, 2026-11-19). With that adapter skipped, WG4's 3D self-test (`app_window_3d`) failed on
  CI: its one case passed the WARP comparison (16 assertions) and then skipped the GPU half, and Catch2 exits
  non-zero when every case it ran skipped. WG4 is now two cases, WARP and GPU, so CI passes the WARP case and
  skips the GPU case. Nothing is loosened: both comparisons are unchanged, and the GPU case still runs on a
  hardware adapter (RTX machine: passed, 16.4 s).
- MG8 full suites (2026-11-19, on 26h's f33de73, MSVC 14.51, RTX 5060 Ti): Release **229 of 229** passed. Debug
  207 of 208: KG6's redraw ratio (enhanced ≤ 1.10× classic) failed at 1.11×, and on a rerun of the test once
  passed and once failed at 1.19×. Owner ruling: the KG6 redraw bound is checked in Release only, like its edges
  bound; Debug records the ratio.
- MG9, hand mutants (as amended by MQ8 and MQ9): **8 of 8 killed**. Each was run alone on a clean rebuild, with its
  target tests passing on unmutated code both before and after.

  | Mutant | Killed by |
  |---|---|
  | Electrodes left matte | the MG2 scene test, the MG3 image test |
  | Glass term power 3 → 1 | the enhanced solid's reference |
  | Glass on opaque insulators | the MG2 scene test |
  | Shadow drawn looking up | the MG5 test |
  | Metal shading in classic | the MG1 bit-identity check, the classic references |
  | Mesh-edge alpha left at 0.75 | the MD6 test |
  | Colourbar change in classic | the classic references with colourbars |
  | Glass rim in dark ink | the edges test, the MG3 glass check |

  A first run of the mutants was discarded. Restoring each file kept its old timestamp, so the build did not
  recompile it and the mutants accumulated. The runner now touches each file it restores and checks a baseline
  before and after.
- **Open:** MG8, the full Release and Debug suites, which wait for the owner's word (2026-11-18). Then the ARCHITECTURE
  entry.
- **Seen, not changed:** in 3D, contours are drawn on top of everything, so they show through the gate and the
  pads. This predates 26g and is the same in classic.

### Unit 26h — scale and orientation (proposed 2026-11-18; on `app/scale-orientation`)

**Owner request (2026-11-18):** as for 26g, items 6 (an axis box with µm ticks) and 7 (a bigger orientation
widget with a click-to-snap view cube, and smooth camera moves). The branch is stacked on 26g's.

**Standing conditions:** as for 26g. These are enhanced-only overlays, except the camera commands, which work in
both looks.

**Today:**
- 2D views have axes with 1-2-5 ticks in µm (`Renderer` draws them with Direct2D from the scene's bounds and
  `ColorScale`-style tick generation).
- 3D views have only the axis triad: three 30 px lines in the lower-left corner (`FieldScene::axis_triad`).
- The camera moves by drag, the keyboard's orbit and zoom commands, and **Fit** (`viewport_fit`, instant).
- There are no named views (26d listed camera presets as not in that unit).

**Design:**

| # | Change | Where | Behaviour |
|---|---|---|---|
| JD1 | **Axis box (item 6).** On enhanced 3D views, the bounds' box is drawn as thin lines (0.55/0.58/0.62, 1 px, depth-tested) on its three faces turned **away** from the camera, ParaView's "grid axes" style, so it frames the device and never covers it. Along three edges (the lower front x edge, the front y edge on the left, the lower z edge on the side) go 1-2-5 ticks in µm, labelled with `Segoe UI` 11.5 px. Each axis is named "x (µm)", "y (µm)", "z (µm)". | Pure `render::axis_box(scene, w, h)` in a new `render/axes3d.hpp`: the edges, the ticks' world positions, their pixels and labels; the renderer draws what it returns | Edges and ticks follow the camera: each axis takes the box edge nearest the viewer on the side away from the device's centre. Labels closer than their own width plus 6 px to a drawn label are dropped. |
| JD2 | **View cube (item 7).** On enhanced 3D views, the triad becomes a 72 px cube in the lower-left corner. It is drawn with Direct2D from the camera's rotation, faces shaded with the 26f lights and labelled Front (+z), Back (−z), Right (+x), Left (−x), Top (+y), Bottom (−y), with x, y, z arrows from its corner. Classic keeps the triad. | `renderer.cpp`; a pure `render::view_cube(scene, w, h)` (the face, edge and corner polygons in pixels) | `nitcad plot` draws it too, so the app and the command give the same image (MG8). |
| JD3 | **Snapping.** In the app, a click on the cube's face snaps to that view (azimuth and elevation from a table: Front 0°/0°, Right 90°/0°, Back 180°/0°, Left 270°/0°, Top 0°/90°, Bottom 0°/−90°). A click on a corner snaps to that isometric view (±45° azimuth, ±35.26° elevation). Zoom and pan reset as in **Fit**. | Pure `render::view_cube_hit(scene, w, h, x, y)` → `std::optional<CameraPreset>`; `main_window.cpp` on a click (not a drag) over it | The cursor becomes a hand over the cube. Clicks elsewhere keep Probe, Measure and Cut line as today. |
| JD4 | **Named views.** **View ▸ Front / Back / Left / Right / Top / Bottom / Isometric**, with the numpad keys as CAD programs use them: 1 front, Ctrl+1 back, 3 right, Ctrl+3 left, 7 top, Ctrl+7 bottom, 0 isometric (the 26f default camera, 35°/25°). Both looks. `nitcad plot` needs nothing new: `--azimuth` and `--elevation` already name any view. | `commands.cpp`, `main_window.cpp` | The equivalent command gives the preset's azimuth and elevation. |
| JD5 | **Smooth camera moves.** Fit, the named views and cube snaps move over 250 ms, eased in and out (smoothstep), taking the shorter way round in azimuth. Elevation, zoom and pan are interpolated linearly in the same ease. A pure `app::camera_between(from, to, t)` gives each frame, and a 16 ms `WM_TIMER` (the playback timer's mechanism) drives it. A drag, a wheel step or a key during a move ends it at the current frame. | `view_settings` (pure), `main_window.cpp` | Instant when Windows' "Show animations" is off (`SPI_GETCLIENTAREAANIMATION`) and in the self-test. The final camera equals the target exactly. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| JQ1 | Where the axis box's ticks go | (a) On the edges that follow the camera (nearest the viewer, on the far side from the centre), ParaView style. (b) Fixed edges at the box's origin corner. | **(a)**: the labels stay readable and outside the device at every angle. Fixed edges end up behind the device or edge-on. |
| JQ2 | The axis box's faces | (a) Only the three faces away from the camera, no gridlines across them. (b) Those faces with gridlines at the ticks. (c) The whole box. | **(a)**: a frame and a scale with no lines across the field. (b) adds a grid behind translucent oxide that reads as structure. (c) covers the device. |
| JQ3 | Is the axis box on by default? | (a) Yes, part of the enhanced look, with an inspector check box **Axis box** (saved, like 26f's) and `nitcad plot --axes on\|off` (default on). (b) Off by default. | **(a)**: scale is information a TCAD view should show unasked. Both switches are needed for screenshots without it. |
| JQ4 | The view cube in `nitcad plot` images | (a) Drawn, as an orientation indicator (the triad's role), so app and command images stay identical. (b) App only. | **(a)**: (b) breaks "same options, same image" (WG4) unless the app draws an overlay outside the renderer. |
| JQ5 | Smooth moves | (a) 250 ms with smoothstep, instant when the OS turns animations off. (b) Always instant. (c) A setting. | **(a)**: it shows where the camera went, and respecting the OS setting covers those who want no motion without a new option. |
| JQ6 | Keys for the named views | (a) The numpad, as in CAD programs (1, 3, 7, 0, with Ctrl for the opposite side). (b) Ctrl+Alt+letters. | **(a)**: the convention users already know. The numpad's digits are unbound in the app today. |

**Owner ruling (2026-11-18):** ruled under the standing instruction: JQ1 (a), JQ2 (a), JQ3 (a), JQ4 (a), JQ5 (a),
JQ6 (a). Implementation waits for the owner's go-ahead, after 26g.

**Acceptance criteria (frozen before code):**

| # | Criterion |
|---|---|
| JG1 | Classic is bit-identical (every existing reference, triad included). The axis box and cube are drawn only in enhanced 3D. |
| JG2 | `axis_box`: for the default camera on a 2 × 1 × 0.8 µm box, the three faces away from the eye are chosen (checked against the face normals' sign along the view direction). Each axis's ticks are 1-2-5 steps covering the bounds, in µm, at least 3 and at most 8 per axis. Each labelled tick lies on its edge, outside the device's projected silhouette. After rotating 90°, 180° and 270° in azimuth and to ±60° elevation, the edges follow the rule of JQ1. No two kept labels overlap. A bounds of zero thickness along z gives no z axis. |
| JG3 | `view_cube` and `view_cube_hit`: at the default camera, the visible faces are Front, Right and Top. A point inside each face's polygon hits its preset. A point near each visible corner hits that corner's isometric view. A point outside the cube gives nothing. The cube's x, y and z arrows match the triad's directions (dot product of the projected directions > 0.99). |
| JG4 | Named views and snaps: each preset gives its azimuth and elevation exactly, with zoom 1 and no pan. The equivalent command names them (`--azimuth`, `--elevation`). The numpad keys and View menu items run them. Top and Bottom draw without a degenerate camera (the world-up switch of `make_view`). |
| JG5 | `camera_between`: t = 0 and 1 give the endpoints exactly; azimuth takes the shorter way (350° → 10° passes 0°, never 180°); smoothstep is monotone. Window tests: Fit and a preset animate (intermediate frames differ from both ends) and land exactly. A drag during a move stops it. With animation off (an override for tests) the move is instant. |
| JG6 | The app: a click on the cube's Top face makes the view Top after the move. A drag that starts on the cube orbits as before. A click elsewhere in Probe mode still probes. The **Axis box** check box defaults on, is saved, is hidden on 2D, and gives `--axes off` in the equivalent command when off. `nitcad plot` with no `--axes` equals `--axes on`. |
| JG7 | Two new WARP references (enhanced): the 3D solid with the axis box and view cube at the default camera, and at a named view (Right). Both pass the existing rules. |
| JG8 | Cost (hardware adapter, KG6 method): drawing the axis box and cube adds under 1 ms. A smooth move's frames each take under 16 ms on `examples/mosfet_3d.deck`. |
| JG9 | Hand mutants, each killed by a test: the box's faces toward the camera chosen instead; labels not de-overlapped; ticks in cm; a cube face hit mapped to the opposite preset; corners ignored; the azimuth's long way round; the endpoint not exact (t = 1 left at the last eased frame); animation ignoring the OS setting; the axis box drawn in classic; `--axes` default off. |

**As built, 26h·1 (2026-11-19).** JD1's "the box edge nearest the viewer" is read per axis. x and z take the nearer
of their two outline edges. y takes the one further left on the screen: JD1's own example, "the front y edge on the
left", is not the nearer one on a box wider than it is deep, such as the test MOSFET. A face seen exactly edge-on
counts as turned away on its lo side. Labels are de-overlapped across all three axes, as JG2's "no two kept labels"
requires. `axis_box` lives in `renderer.cpp` beside `make_view`; it is pure, with no GPU work.

The three 26g enhanced references now change by design, because the axis box is on by default in enhanced (JQ3). As
in 26g, they are rebuilt in 26h·4 only after the owner sees the images.

**As built, 26h·2 (2026-11-19).** `view_cube` and `view_cube_hit` (axes3d.hpp) are pure, like `axis_box`. The cube
is 72 px across its diagonal, centred 88 px in from the view's lower-left corner, and always orthographic. A face's
name is drawn only where the face is wide enough for it (6 px a character) and at least 16 px tall. At the default
camera only "Front" fits; "Top" and "Right" are too foreshortened. Corners come before faces in hit testing, within
7 px. A corner snaps to its isometric view: azimuth atan2(x, z), elevation ±35.26°. The plain Isometric view is
the default camera (35°, 25°). A snap resets zoom and pan as Fit does. It is instant for now; 26h·3 makes it smooth.
A cube click takes priority over Probe, Measure and Cut line. A drag that starts on the cube orbits as anywhere else.
The cursor becomes a hand over the cube.

**As built, 26h·3 (2026-11-19).**
- **Where the named views sit.** They are in a **View ▸ Loo&k from** submenu (&Front, &Back, &Left, &Right, &Top,
  B&ottom, &Isometric), not as seven flat View items. WG5 requires a unique mnemonic for every item, and View had
  only j, k and z free. The commands test now checks mnemonics within each submenu, and the submenu's own mnemonic
  within its menu; the rule itself is unchanged.
- **Where the numpad keys live.** The view's own key handling takes them (`viewport.cpp`), as the arrows, + and −
  are, not the accelerator table. The menu shows them ("Num 1", "Ctrl+Num 1", ...), and the deck editor still types
  numpad digits.
- **Who animates.** `WindowOptions::animate` is off by default and on only in the application (`app/main.cpp`), so
  the tests and the self-test land at once. `set_animate(true)` forces moves on whatever Windows says, so the
  animated window test also runs on a machine with animations off.
- **What ends a move.** Orbit, pan and zoom (drag, wheel, key) end a move at its current frame. So does a change of
  tab, which stops it at its next tick.
- **Fit.** It animates zoom and pan only; the angles stay.
- **JG8.** The slowest frame of a smooth move drew in 13.8 ms on WARP (69 samples). The hardware figure is for
  26h·4.

**Results, 26h·4 (2026-11-19).**
- **JG1.** Every classic reference passes unchanged; only the three 26g enhanced references differ (1.1% of pixels,
  the box and cube). Rebuilt after the owner approved the images (2026-11-19), byte-identical to the ones shown.
- **JG7.** `mosfet3d_solid_right_enhanced` is new, a VG4 section. Its render showed the cube's arrows ending inside
  the cube on head-on views. The arrows now reach 10 px past the cube, each label 8 px past its arrow's tip, and the
  cube sits 88 px in from the corner.
- **JG8**, by KG6's method on `examples/mosfet_3d.deck`, medians of 50:

  | Adapter | Box and cube cost (bound 1 ms) | Slowest move frame (bound 16 ms) |
  |---|---|---|
  | NVIDIA GeForce RTX 5060 Ti | +0.13 ms (3.60 vs 3.47 ms) | 5.07 ms |
  | WARP, recorded only | 15.99 vs 14.26 ms | 19.2 ms |

  JG8 is met.
- **JG9.** All 10 hand mutants were killed by tests. The tree was restored and touched after each mutant, with
  baselines passing before and after.

  | Mutant | Killed by |
  |---|---|
  | faces toward the camera | `render axes3d` |
  | labels not de-overlapped | `render axes3d` |
  | ticks in cm | `render axes3d` |
  | face hit to the opposite preset | `render axes3d` (cube) |
  | corners ignored | `render axes3d` (cube) |
  | the azimuth's long way round | `app views 26h` (camera_between) |
  | endpoint not exact | `app window 26h: smooth` |
  | OS setting ignored | `app views 26h` (smooth_moves) |
  | box drawn in classic | VG4 "3D solid" |
  | `--axes` default off | `app views 26h` (plot images) |

  The OS-setting rule was moved into the pure `app::smooth_moves`, so a test can reach it. The first try of that
  mutant failed only on an unused-parameter warning and was rerun with the parameter used.
- **Still open.** The full Release and Debug suites wait for the owner's word.

**Slices (each one local commit):** 26h·1 `axis_box`, drawn, its switch (JD1; JG2, JG6 part); 26h·2 `view_cube`, hit
testing, snapping (JD2, JD3; JG3, JG6 part); 26h·3 named views and smooth moves (JD4, JD5; JG4, JG5, JG8); 26h·4
references, suites, mutants (JG1, JG7, JG9), the ARCHITECTURE entry.

---

### Unit 61 — restart of transient, trace and small-signal runs; `initial`; time and frequency suffixes (proposed 2026-11-19; on `io/restart-transient-trace`)

**Owner request (2026-11-19):** "start Restart of transient and trace runs, starting a run from a checkpoint
(initial), unit suffixes on times", the 11.3 row "Input/output, checkpoint/restart", next step "Restart of transient
and trace runs; a run started from a checkpoint (`initial`); unit suffixes on times" (also "Not in Unit 43", 6.5).
The branch starts from `main` `8f84d04`. Unit 60 (PR #96) is not in it; the two touch different code.

**Starting point (read from the code):**
- **Sweep restart (Unit 31, DK9):** `[sweep] restart = <checkpoint>` finds the plan point whose bias the checkpoint
  holds, skips the points up to it and starts `sweep_bias` from the checkpoint's fields. The run identity is not
  checked. The restarted run reaches the uninterrupted points to 1e-9 relative, not bit for bit, because the state
  passes through physical units.
- **A checkpoint** (`io::Checkpoint`, content kind 7) holds the mesh digest, run identity, bias, thermal bias, node
  fields and provenance, and nothing about the analysis that wrote it.
- **Transient:** `checkpoint =` with `[transient]` is refused (AQ4: "its end state is not a steady point"). The
  integrator carries up to three accepted states (scaled x, stored charge, trap occupancy, contact charge, time),
  the next step size, and restarts its history at waveform breakpoints (`solve/transient.cpp`). A restart from
  fields alone would start over with backward Euler and a new step size, so it would not continue the same run.
- **Trace:** the continuation carries the scaled state with λ, the unit tangent (oriented by the previous
  tangent's arc row), the arc step `ds` and the bordered row's entry `k`. Past a fold the bias moves away from
  `end_V`, so a restart that re-derives the direction from `end_V` would turn back onto the branch already traced.
- **`initial`:** `solve_transient`, `trace_bias` and `solve_small_signal` already take an initial state (a starting
  guess). Only the deck lacks a way to give one.
- **Times** are bare seconds (DK5 kept, AQ2) and frequencies bare hertz. Lengths take `cm`, `um`, `nm`.

**Standing conditions:** opt-in. Without `restart` or `initial` every analysis is bit-identical to `main`, and so is
every file written by an analysis that existed before (RG1). Old checkpoints load unchanged (IO2). Bare numbers keep
their meaning. Criteria are fixed before code. Release and Debug suites with `/W4 /WX`, plus hand mutants. Built
with MSVC 14.51. Nothing is pushed or merged without authorization.

**Left out of this unit:** extending a finished run (a longer `t_end`, a further `end`), because it changes the step
defaults derived from `t_end`; restart under `[adapt]` (refused as today); periodic checkpoints during a run
(RQ5); voltage suffixes (`mV`); the GUI's restart controls.

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| RD1 | 1 | **Suffixes.** Times take `s`, `ms`, `us`, `ns`, `ps`, `fs`; frequencies take `Hz`, `kHz`, `MHz`, `GHz`, `THz`. They are read by shifting the decimal exponent of the text, so `2.5ns` parses to the same double as `2.5e-9`. They apply where the option or argument is a time or a frequency: `t_end`, `dt_*`, `output_times`, the time arguments of `step`, `ramp`, `pulse`, `pwl` and `sine` (`f` a frequency), and `[small_signal] frequencies`, `logspace` included. A suffix of the wrong kind (`1ns` as a voltage, `1MHz` as a time) is an error with its line. | Bare numbers unchanged (RG1). |
| RD2 | 2 | **`initial = <checkpoint>`** in `[transient]`, `[small_signal]`, `[trace]` and `[sweep]`: the checkpoint's fields become the run's starting guess (the existing `initial` argument). It is refused with `restart` in the same section, on another mesh (as `load_checkpoint`), and for a checkpoint holding resume state (RD3; use `restart`). | A guess only: the bias need not match (RQ3). |
| RD3 | 3 | **Transient resume state.** `results::TransientResume`: the integrator's accepted history (up to three states as scaled vectors, bit for bit, with their times), the next step, the number of accepted steps, the next output-time index and the breakpoint state. `solve_transient` returns it at the end of every run and takes an optional one to continue from. Continuing reproduces the uninterrupted run's arithmetic. | Absent: today's code path (RG1). |
| RD4 | 3 | **Checkpoints with resume state.** `io::Checkpoint` gains an analysis kind (absent = steady; IO2) and an optional resume record. `[transient] checkpoint =` is allowed and writes the last accepted step with its resume state (AQ4 revised). `[transient] restart =` checks the run identity: **equal identity required** (RQ2), because the state belongs to that run. The restarted run's files hold the steps after the checkpoint, its first point the checkpoint's step. | |
| RD5 | 4 | **Trace resume state.** `results::TraceResume`: the scaled state with λ, the tangent, `ds`, `k` and the point count. `trace_bias` returns and accepts it the same way. `[trace] restart =` and `checkpoint =` work as RD4. | |
| RD6 | 5 | **Small-signal restart.** `[small_signal] restart = <checkpoint>` continues after the operating point the checkpoint holds, as the sweep's (DK9), at every frequency of the remaining points. | |
| RD7 | 6 | **Records:** ARCHITECTURE 6.5 and 6.6 "As built", 11.3; DECISIONS results. | |

**Criteria (frozen once approved):**

| # | Gate |
|---|---|
| RG1 | Bit-identity: every existing test passes unchanged. Every example and `decks/` deck that runs in under 60 s writes a CSV table byte-identical to `main`'s, and every test that compares a deck run with the C++ call bit for bit still passes. |
| RG2 | Suffixes: every suffix gives the same double as the bare exponent form (`2.5ns` ≡ `2.5e-9`, `10GHz` ≡ `1e10`), for each key and each waveform argument; a deck written with suffixes gives the bare deck's run identity; the wrong kind and an unknown suffix name the line. |
| RG3 | `initial`: a transient, small-signal, trace and sweep run from a checkpoint of its own starting point equals the run without it to 1e-10 relative, with no more Newton iterations at the first point. Refused: with `restart`, on another mesh, and with a resume checkpoint. |
| RG4 | Transient restart: on the hand-built diode with a pulse (BDF2, adaptive), runs stopped after 1, 2 and 3 accepted steps, after a breakpoint and mid-ramp each restart, and **the joined run equals the uninterrupted run bit for bit**: times, steps, orders, error ratios, currents and snapshots. The same holds with interface traps, and with electrothermal on. A restart with another identity, a steady checkpoint, or the run's final step is refused. |
| RG5 | Trace restart: a trace stopped after k points (k = 1 and a k past a fold, on Unit 19's fold case) and restarted equals the uninterrupted trace bit for bit from point k on. A steady checkpoint and another identity are refused. |
| RG6 | Small-signal restart: the remaining operating points equal the uninterrupted run's to 1e-9 relative (as DK9), at every frequency; a checkpoint at the last point is refused. |
| RG7 | Files: Unit 29's and Unit 42's checkpoint fixtures load unchanged as steady checkpoints; a resume checkpoint round-trips bit for bit; a damaged resume record is refused. |
| RG8 | Deck and CLI: every new mistake names its line; `nitcad run` on a stopped transient (exit 2) writes a checkpoint that a restart completes (exit 0), and likewise for a trace. |
| RG9 | `/W4 /WX`. Hand mutants (one history state dropped, the next step not restored, the breakpoint state lost, the tangent's sign flipped, `ds` not restored, the identity check removed, the suffix exponent off by 3, the wrong suffix kind accepted, a skipped small-signal point), each caught. |

**Size (estimate):** L, about 1.5–2k lines, six slices. Slices 3 and 4 touch the transient and trace solvers.

**Questions for the owner (recommended option first):**

| # | Question | Options |
|---|---|---|
| RQ1 | Restart exactness | **Bit-identical continuation through solver resume state (recommended).** / Restart from fields only, to a tolerance (simpler, but a transient would restart its integrator and change every later step). |
| RQ2 | Run identity on restart | **Transient and trace: equal identity required. Sweep and small signal: as DK9, bias match only (recommended).** / Required everywhere (changes Unit 31's sweep restart). |
| RQ3 | `initial` with a bias other than the run's start | **Allowed; it is only a guess (recommended).** / Refused. |
| RQ4 | Which analyses get `initial` | **Transient, small signal, trace and sweep (recommended).** / Transient and small signal only (the 11.3 wording). |
| RQ5 | When a transient or trace checkpoint is written | **At the end of the run, completed or stopped, as sweeps today (recommended).** / Also every N steps (`checkpoint_every`), for crashes. |
| RQ6 | Suffix scope | **Times and frequencies (recommended).** / Times only. |
| RQ7 | Small-signal restart | **In this unit (recommended).** / Later. |
| RQ8 | The restarted run's output files | **Hold the steps after the checkpoint, the first being the checkpoint's own step (recommended).** / Joined to the earlier file. |

**Owner ruling (2026-11-19, standing instruction "always choose recommended", 2026-11-11):** RQ1–RQ8 take the
recommended option. Pushing, the PR and the full Release/Debug suites still wait for the owner's word.

**Results (2026-11-19, Release, MSVC 14.51; the full suites not yet run, on the owner's word):**

| Gate | Result |
|---|---|
| RG1 | The 5 tracked `decks/` and `examples/` decks that run in under 60 s write CSV tables byte-identical to the branch base `9225e6f`'s (built in a separate worktree). The 5 FinFET decks take over 60 s with both builds and are outside the gate. The deck–C++ bit-for-bit tests pass. Existing expectations changed only where this unit changes behaviour by ruling: `[transient] checkpoint =` is accepted (AQ4 revised, RD4), and `restart` is accepted in `[transient]`, `[trace]` and `[small_signal]` (`tests/deck/analyses_test.cpp`). |
| RG2 | Units: 3 cases pass. `2.5ns` gives the same double as `2.5e-9`, and the wrong kind, an unknown unit and a unit on a voltage each name their line. |
| RG3 | `initial`: 2 cases pass. The first state takes 6 → 1 Newton iterations in all four analyses. |
| RG4 | Transient resume is bit for bit. The diode pulse (BDF2 adaptive) resumes at k = 1, 2, 3, at the breakpoint, after it, mid-ramp, at the top of the rise and at the next-to-last point. The trapped MOS capacitor and the electrothermal diode also resume bit for bit. Refusals: another run, the final step, `initial` with resume, nonlocal tunnelling. Deck: stopped after 3, 20 and n/2 steps, then restarted, equals the uninterrupted run. |
| RG5 | Trace resume is bit for bit. The diode resumes at converged attempts 0, 1, 2, n/2 and n−1. The electrothermal runaway fold resumes before, at and past its top. Stops on rejected attempts resume from the last accepted point. Refusals: another run, a finished trace, an initial state, a state that does not fit. Deck: restarted after attempts 1 and 5 equals the uninterrupted trace. A steady checkpoint and another trace are refused. |
| RG6 | Small-signal restart from a sweep's checkpoint and from a small-signal checkpoint: the 2 remaining points' operating currents and admittances, at 1 MHz and 1 GHz, equal the uninterrupted run's to 1e-9 relative. The last point, a bias that is not a point and a transient checkpoint are refused. |
| RG7 | Checkpoints: the transient and trace resume records round-trip bit for bit. The pre-Unit 42 fixture loads as steady with no resume state. A damaged file is refused. |
| RG8 | `nitcad run` on a transient and on a trace, stopped by Ctrl-Break (exit 2), writes a checkpoint, and a restart completes it (exit 0); this also passes under ctest. A run that stops on its own solver limits (for example `max_iterations = 1`) is correctly refused as a restart with other limits, because the limits are part of the identity. New deck mistakes name their line. |
| RG9 | `/W4 /WX` clean. 10 of 10 hand mutants caught: a history state dropped (2959 failures), the next step not restored (1856), the landing state lost (4), the tangent's sign flipped (2), `ds` not restored (1), the transient identity check removed (1), the trace identity check removed (1), the `ns` exponent off by 3 (10), a frequency unit accepted as a time (1), a skipped small-signal point (1). |

Also fixed: a stopped small-signal restart's message counted the skipped points among those left to solve.

---

### Unit 26i — fast camera redraw (proposed 2026-11-19; on `app/redraw-speed`)

**Owner request (2026-11-19):** rotating and zooming `decks/finfet_3d_shaped_swept.deck` in the app lags; make
it faster "without decreasing its quality". The branch starts from `main` `344dd79`.

**Today (read from the code):**
- **Every camera move rebuilds the frame from nothing.** A drag's `WM_MOUSEMOVE` (`viewport.cpp`) calls
  `viewport_orbit` → `present_camera` → `Viewport::present` → `Renderer::present` (`renderer.cpp`), which
  (1) waits for the GPU to go idle, (2) in `render_field` copies every triangle batch into a new vertex array
  and turns every line into a two-triangle quad, (3) creates new upload buffers and new default-heap vertex
  buffers (`CreateCommittedResource` per frame) and copies everything to the GPU, then (4) draws, and
  `Present(1, 0)` waits for the vertical blank. Only the camera changed: `apply_camera` writes azimuth,
  elevation, projection, zoom and pan, and nothing else in the scene.
- **The comment says otherwise:** "the vertices live in default-heap buffers … (a window keeps them across
  frames)", but nothing keeps them.
- **One present per mouse message.** A mouse reports at 125–1000 Hz and each present blocks on the vertical
  blank, so the messages queue and the view trails the cursor.
- **`update_metrics` per frame** scans every mesh node for the bounds and the field for its largest value, and in
  the quality views recomputes every element's angles (`render::element_quality`).
- `Renderer::frame_time` (26f's KG6) already measures a redraw on resident geometry. It re-executes the recorded
  passes, so it times the GPU, not the CPU rebuild above.

**Standing conditions:** no change to what is drawn. `Renderer::draw` (images, `nitcad plot`, every WARP
reference) is untouched, so every reference passes by construction. The window's frames equal what a full
rebuild draws, **bit for bit** (ND1, NG1). The 2× supersampling, the transparency passes, the feature edges, the
finishes, the ground shadow and the axis box stay as they are. Criteria are fixed before code. Release and Debug
suites with `/W4 /WX`, plus hand mutants. Built with MSVC 14.51. Nothing is pushed or merged without
authorization.

**Left out of this unit:** overlapping CPU and GPU work across frames (two frames in flight, NQ3); a lower level
of detail or quality while dragging (NQ2: refused, it lowers quality); faster scene building after a field or
option change (the wait after "Building the deck" or a field change is another path).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| ND1 | 1 | **Resident geometry.** The window target keeps the GPU vertex and line buffers (and the draw ranges) of the last full frame. `Renderer::present_camera(scene, target, w, h)` redraws with them, building only the per-frame constants (the camera matrices, the lights' eye-space vectors, the shadow plane) and the camera-dependent axis box lines (JD1 picks the faces turned away from the camera; they get a small buffer of their own). `Renderer::present` is unchanged and refreshes the resident buffers. | The constants come from the same code as today, so the bits match. |
| ND2 | 1 | **Who calls it.** `MainWindow::present_camera` is the only caller: the orbit, pan, zoom, keyboard, Fit, named-view and smooth-move paths, which change nothing but the camera (NQ1). Every other present (a new field, look, clip, cut, edges, contours, colormap, a resize) takes the full path. A target with no resident geometry, or with another scene's, falls back to the full path. | The scene gains a revision counter that every edit outside `apply_camera` bumps. The target stores it, and a mismatch means the full path. |
| ND3 | 1 | **Debug check.** In Debug builds `present_camera` hashes the scene's vertices and lines and compares the hash with the one stored at the full frame. A mismatch is a failed assertion (a missed revision bump). | Release skips the hash. |
| ND4 | 2 | **Coalesced input.** A mouse move, wheel notch or camera key updates the camera and calls `InvalidateRect`. The draw happens in `WM_PAINT`, which Windows sends only once the message queue is empty, so a burst of moves gives one frame at the latest camera. The smooth moves (JD5) keep their timer. | No camera step is lost: the camera takes every message. |
| ND5 | 3 | **Metrics cached.** The bounds and the node and element counts are computed once per workspace, the largest value once per field, and the quality summary once per quality scene. `update_metrics` then only formats the text and the "Drawn in" time. | The text is unchanged. |
| ND6 | 3 | **Records:** ARCHITECTURE 6.14 "As built" (the window's resident path); DECISIONS results with the before and after times. | |

**Criteria (frozen once approved):**

| # | Gate |
|---|---|
| NG1 | **Bit-identity.** On WARP, for a 3D scene in each look (classic; enhanced with metal and glass; clipped; exploded; with edges, contours and the axis box), a `present_camera` frame after a camera change reads back **byte-identical** to a full `present` at the same camera, over 8 camera poses including both projections and an elevation below −5° (no shadow). Every existing reference test passes unchanged. |
| NG2 | **Invalidation.** Each scene edit (field, scale, colormap, look, finish, clip, cut, exploded, edges, contours, axis box, transparency, resize) followed by a camera move draws the new scene, byte-identical to a full present. A target holding another scene's buffers takes the full path. |
| NG3 | **Coalescing.** In a window test, 200 drag `WM_MOUSEMOVE` messages posted at once give at most 3 presents, and the final camera equals the one that 200 separate moves give. The JD5 moves, Fit and the named views still end on their exact poses. |
| NG4 | **Speed** on `decks/finfet_3d_shaped_swept.deck` and `examples/mosfet_3d.deck`, enhanced look, the hardware adapter, 1280×960: the time of a camera-only present (`last_present_ms`), median of 50, is **at most a third** of `main`'s before this unit (both measured in this unit and recorded). WARP times are recorded only. |
| NG5 | **Memory.** The resident buffers are freed with the target and replaced, not added to, on each full present: the GPU memory after 100 field changes equals that after one, within 1%. |
| NG6 | **Metrics:** the status text equals `main`'s for every view kind, the quality views included. |
| NG7 | `/W4 /WX`. Hand mutants (the revision not bumped on a clip change, the constants not rebuilt, the axis box not rebuilt, the resident buffers kept after a resize, the `WM_PAINT` coalescing bypassed, the cached bounds kept across workspaces), each caught. |

**Questions for the owner (recommended option first):**

| # | Question | Options |
|---|---|---|
| NQ1 | How the renderer knows only the camera moved | **An explicit camera-only call from the one window path, guarded by a scene revision and a Debug hash (recommended).** / Hash the whole scene every frame (safe without the app's help, but hashing tens of MB per frame costs much of the saving). |
| NQ2 | A cheaper frame while dragging (no supersampling, no edges) | **No: the frame stays identical (recommended).** / Yes, with full quality when the drag stops. |
| NQ3 | Two frames in flight (the CPU builds the next while the GPU draws) | **Later, if NG4's numbers call for it (recommended).** / In this unit. |
| NQ4 | Where the input is coalesced | **`WM_PAINT` via `InvalidateRect` (recommended).** / A 16 ms timer (adds up to a frame of latency). |

**Owner ruling (2026-11-19, standing instruction "always choose recommended", 2026-11-11):** NQ1–NQ4 take the
recommended option. Pushing, the PR and the full Release/Debug suites still wait for the owner's word.

**Size (estimate):** M, about 400–600 lines with tests, three slices. Slice 1 touches `renderer.cpp` and
`viewport.cpp`; slices 2 and 3 the app.

**Results (2026-11-19, Release, MSVC 14.51; `690d8bb` slice 1, `bd0364e` slices 2–3; Debug only for NG1, then
stopped on the owner's word "dont start debug until i say"; the full suites not run):**

| Gate | Result |
|---|---|
| NG1 | **Met.** 5 scenes (classic solid; enhanced solid with metal, glass, shadow and axis box; enhanced potential with mesh edges and contours; enhanced clip; enhanced exploded) × 8 poses (both projections, elevations −30° to 89°, zoom and pan): every `present_camera` frame, read back with `Renderer::last_image`, is byte-identical to `draw` at the same camera, and each reused the kept buffers. Release and Debug (316 checks; Debug run before the owner's word). All 98 render test cases pass unchanged: every reference matches. |
| NG2 | **Met.** Colormap, classic and enhanced look, a finish, transparency, axis box off, shadow off and new field values, each followed by a full present and a camera move: the new scene, byte-identical, reused. Another scene object and the same object with another batch count take the full path. A resize reuses the buffers, and a plot in between drops them. |
| NG3 | **Met.** In a shown, input-disabled window, 200 posted drag moves give 1–3 frames (5 runs), at the camera that 200 moves drawn one by one give, to the bit; the frame is `draw`'s. 20 posted wheel notches: at most 3 frames, every notch taken. The window had to be input-disabled: with the cursor over it, the desktop mouse added moves during the drag's capture. |
| NG4 | **Met by the owner's ruling (2026-11-19, "choose recommended"): judged on the deck that lagged, `mosfet_3d` recorded.** RTX 5060 Ti, 781×616, enhanced, medians of 50 (`--orbit-benchmark 50`, camera-only then full on each frame): `finfet_3d_shaped_swept` **6.72 ms against 92.72 ms (13.8×)**, met; `mosfet_3d` **2.25 ms against 3.81 ms (1.7×)**, short of a third. Timing the stages of a camera-only `mosfet_3d` frame: building the constants 0.33 ms, the passes 0.18 ms, Direct2D 0.27 ms, the copy and `Present(1, 0)` 1.34 ms. The present is paced by the vertical blank in a loop that draws two frames per step, so the gate's ratio cannot fall to a third for a scene this small. |
| NG5 | **Met.** WARP, the process's private bytes after 100 field changes, each presented in full and then turned: 0.988–1.000 of those after the first (3 runs). WARP keeps freed memory: the first hundred draws grow the process by some 100 MB whatever is kept, so the test warms up to that first. At 300 changes the bytes were below those at 100: nothing accumulates. |
| NG6 | **Met.** The strip on the 3D MOSFET: unchanged over camera moves for the potential, the electrons and mesh quality; the electrons' differs from the potential's; back from quality to the field, the field's again; on opening the 2D MOSFET, its own size, triangles, no tetrahedra, and XG7's quality figures. After each field change the window's frame is the scene's, drawn in full. |
| NG7 | **Met.** `/W4 /WX` clean. 8 of 8 hand mutants caught: `Viewport::show` taking the camera path (2 failures), the batch sizes not compared (crash), the axis box not rebuilt (crash), the kept buffers not dropped by a plot (1), the `WM_PAINT` coalescing bypassed (4), no workspace generation bump (4), no field generation bump (1), the camera constants kept from the first frame (2). Two first forms survived as equivalent mutants (one bump removed of several that fire on the same path), and the first form of the coalescing mutant did not compile under `/WX` (unreachable code). They were rerun in the forms above. |

- **NG5 measured on the GPU, on the owner's word "fix it" (2026-11-20).** The private-bytes check was flaky: after
  the hundred-change warm-up (`3e16edc`) it still failed PR 105's CI run (1.013) and one local run in three (1.017),
  because WARP's and the CRT's heaps drift 1–2.5 MB in 100 changes whatever the renderer keeps. NG5 is about GPU
  memory, so the test now measures that: the process's D3D12 memory on the WARP adapter (DXGI
  `QueryVideoMemoryInfo`, local plus non-local), after one change and after a hundred more, within NG5's 1% (not
  loosened). It is the same to the byte (19,202,048 bytes; ratio 1.0 in 5 of 5 runs), the warm-up draws are gone,
  and a mutant that keeps every replaced `ResidentGeometry` alive gives 10.9. `nitcad_render_test` links `dxgi`.
- **Design change, ND2:** no scene revision counter. In the app the scene changes only through `Viewport::show`,
  which presents in full, so the target's scene address and batch sizes, with the Debug hash, are enough. NG2
  tests the cases.
- **Fixed on the owner's word "fix it" (2026-11-19), two faults that were there before this unit:**
  - **Closing a shown `MainWindow` crashed in the test executable** (`3dd312d`). The crash was in releasing the
    taskbar list. The application initializes COM in `main`; the test executable did not, so the list came from
    another thread's implicit apartment, gone by then. `MainWindow` now initializes COM (STA) for its lifetime,
    balanced after that release; `Viewport::release` also stops a swap chain being made again while the window is
    destroyed. A regression test closes a shown window with a deck built; it crashed before the fix. NG3 still
    drives a `Viewport` directly, which is where coalescing lives.
  - **KG6 wandered** (`a148e79`): 1.04–1.18 on one build, because all classic frames were timed and then all
    enhanced ones. The looks now alternate. Interleaved, the test caught slice 1's per-frame axis-box buffer (1.12–1.16
    against 1.06–1.11 before the unit): the renderer now keeps one, and keeps the overlays' DirectWrite layouts
    across frames. That left 1.07–1.10, still at the bound. Profiled: the GPU passes cost the same in both looks
    (about 0.2 ms); the enhanced extra is CPU work, Direct2D for 26g's colour bar and 26h's view cube and axis box
    (about 0.22 ms) and the feature edges' vertices (about 0.07 ms), on a 5 ms draw. **Owner ruling (2026-11-19):**
    KG6's bound is **1.15×**.

---

### Unit 26k — interactive device geometry editor (proposed 2026-11-20; on `app/geometry-editor`)

**Owner request (2026-11-20):** "I want to implement a professional, fully interactive Device Geometry Editor in
NiTCAD. Currently, users cannot create semiconductor device structures directly using the mouse. I want
functionality similar to commercial TCAD software, where users can draw, select, move, resize, and edit device
geometries directly in the viewport." The branch starts from `app/console-log` `a25e0da` (Unit 26j and look A,
unmerged), since the editor sits in the same window code.

**Where things stand.** A device is typed: a tensor deck's regions are a background region plus `box`es in deck
order, each replacing what it covers; a `generate2d` deck's regions are `polygon`s that must *tile* the domain
exactly (shared edges with bit-identical vertices, no overlap, no nesting, the interface rules M7/U12). Contacts are
boundary `patch`es (generated) or `box`es on a named face (tensor). Nothing can be drawn; the Device menu only opens
the deck at its sections.

**The model the editor uses (the commercial one).** Shapes are layered: a background region fills the domain, and
every shape after it replaces what it covers ("replace" mode). That is already the tensor deck's rule. A tiling is
then computed from the layers, never drawn by hand: the user never has to make two polygons share a vertex exactly.

**Standing conditions:** the deck stays the only source of truth (VQ13): the editor reads the deck and writes deck
text, nothing else, through the deck editor's undo; every deck that builds today builds to the same bits (the new
deck key is opt-in); Win32 + Direct2D only, no new dependency; the geometry model has no Win32 type (testable);
criteria fixed before code; targeted Release tests per slice, full suites only on the owner's word; no screen
automation; nothing pushed or merged without authorization.

**Left out (later units):** 3D editing (extruding the cross-section into a swept or generated 3D deck); drawing
doping (Gaussian windows, implants); process emulation (etch, deposit); curved shapes (arcs, as polylines only);
importing GDS/DXF; meshing controls drawn on the canvas (refinement boxes).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| GD1 | 0 | **Layered generated decks.** `[mesh] type = generate2d` gains `layered = on` (default off: today's tiling rule unchanged). Layered, the domain is `[mesh] domain = x0 x1 y0 y1`; one region without shapes is the background; regions give `box = x0 x1 y0 y1` (as tensor decks) or `polygon`, applied in deck order, each replacing what it covers. | Opt-in: no existing deck changes meaning. |
| GD2 | 0 | **The partition kernel** (`mesh::layered_partition2d`, mesh layer, exact): the layers resolved into the tiling `generate2d` takes. Rectilinear shapes in slice 0 (an overlay of every shape's x and y breakpoints, each cell to its topmost shape, cells merged per region and connected component into simple polygons, T-junction vertices inserted on shared edges so every edge is conforming). Slanted edges in slice 4. Holes (a region enclosing another) are cut along a breakpoint line into simple polygons, as `generate2d` takes no nesting. | Coordinates are the inputs' own doubles, never recomputed: shared vertices are bit-identical by construction. |
| GD3 | 1 | **The geometry model** (`app::Geometry`, no Win32): the domain, the shapes (name, material, a rectangle or polygon, order), the contacts (name, kind, segment on the boundary), read from a deck (tensor, generate2d tiled or layered) and written back as deck text by minimal edits (only the `[mesh]` geometry keys and the regions' `box`/`polygon` lines change; comments, order and every other section kept). A tiled `generate2d` deck reads as one shape per polygon (exact; editing it writes it layered). | |
| GD4 | 1 | **Edits and undo:** add, delete, move, resize, set a vertex, insert/remove a vertex, rename, set material, reorder (bring forward, send back), add/move/resize a contact, duplicate, mirror; each an operation on the model, applied as one deck edit (one undo step in the deck editor; Ctrl+Z in the canvas undoes it). | |
| GD5 | 1 | **Snapping and hit testing:** a grid (1, 2, 5 × 10ⁿ µm steps chosen from the zoom), vertices, edge midpoints, edges, the domain boundary; Shift constrains to horizontal/vertical; coordinates quantized to 0.1 nm so decimal deck text round-trips exactly. Hit tests in screen pixels (handles 7 px at 96 dpi, scaled). | |
| GD6 | 2 | **The Geometry tab** (centre, after Start): a Direct2D canvas with rulers in µm, the grid, the shapes filled by material colour (the tree's swatches) in layer order, outlines, contacts as thick coloured segments with their names, the selection with 8 resize handles, rubber-band selection, zoom at the cursor (wheel), pan (middle drag, Space+drag), Fit. The cursor's position in µm in the status bar. | 2D decks only; a 3D deck shows "3D geometry editing is a later unit". |
| GD7 | 2 | **Tools** (a strip on the canvas and the Edit menu, each with a key): Select (V), Rectangle (R), Polygon (P: click vertices, double-click or Enter closes), Contact (C: drag along the boundary), Pan (H), Measure (M, reuses the probe's ruler). Arrow keys nudge by a grid step (Shift: ×10); Delete; Ctrl+C/V/D; Escape cancels. | |
| GD8 | 2 | **Inspector, Geometry section:** the selection's name, material (drop-down of the materials the deck build knows), layer order, x/y/width/height (rectangle) or the vertex list (polygon), contact kind and its keys; typed values apply on Enter. The tree's regions and contacts select on the canvas and back. | |
| GD9 | 3 | **Live validation:** the layout checked as it changes with the build's own checks (the partition, `generate2d`'s interface rules, a contact off the boundary or overlapping another), each problem a marker on the canvas and a line in the console's Messages; Run disabled with the reason while invalid. **Mesh preview** button: builds the mesh in the background (the existing preview path) and overlays it. | |
| GD10 | 4 | **Slanted edges** in the partition (trapezoids: sidewalls, STI): segment intersections in exact arithmetic (the mesh layer's `exact_integer` / orient2d), intersection points rounded once and shared. | |
| GD11 | — | **Commits:** one branch, the proposal, then a commit per slice, then ARCHITECTURE (the app's chapter and 11.x row) and the records. The owner may stop after any slice. | |

**Criteria (fixed before code; a target missed is recorded as a limitation, never loosened):**

| # | Criterion |
|---|---|
| GE1 | **Nothing else moves:** every example and test deck builds to the same bits as before (the existing deck and mesh tests unchanged and passing). |
| GE2 | **Partition exact:** for 1000 seeded random rectilinear layouts (up to 30 shapes, holes included) the kernel's output passes `generate2d`'s tiling validation, the region of every cell centre equals the topmost covering shape's, and the areas per region equal the exact layered areas. A layered deck and its hand-tiled equivalent build to the same mesh bits. |
| GE3 | **Round trip:** every 2D example and fixture deck read into the model and written back with no edit gives byte-identical text; after an edit, only the geometry lines differ (checked by a line diff), and the deck builds. |
| GE4 | **Undo:** a seeded sequence of 200 random edits, each undone in reverse, restores the deck text byte-identical at every step; redo replays it. |
| GE5 | **Speed:** a drag (move or resize) redraws in ≤ 16 ms per frame median on WARP with 200 shapes and 5000 vertices; the partition of 200 rectangles takes ≤ 5 ms; the deck written back per edit ≤ 2 ms. |
| GE6 | **Hit testing:** handles, vertices, edges and contacts picked correctly at 96, 144 and 192 dpi and at zooms 1e-2 to 1e4 px/µm (table test). |
| GE7 | **Validation agrees with the build:** for a table of invalid layouts the editor's message equals the deck build's error, and every valid one builds. |
| GE8 | **Interaction:** window tests drive the canvas by messages (no screen automation): draw a rectangle, move and resize it, draw a polygon and edit a vertex, draw a contact, undo/redo, and the deck text and device are as expected; a drawn MOSFET (substrate, oxide, gate, source, drain, body) runs its sweep. |
| GE9 | **Mutants:** at least 8 per slice caught (e.g. the replace order reversed; a T-junction vertex not inserted; a hole not cut; the quantum dropped; snap to the wrong candidate; undo of a reorder not restoring the order). |
| GE10 | **Layering:** the model in `app` without Win32 headers; the kernel in `mesh` without `app` or `deck`; no new dependency. |

**Questions (each with a recommendation; ruled by the standing "choose recommended" rule, 2026-11-11):**

- **GQ1 — what the deck stores.** (a) the layered shapes (`layered = on`), the user's intent kept; (b) the computed
  tiling as plain `polygon`s, no new key, but a reopened deck has lost which shape was drawn over which.
  **Recommended (a)**, ruled (a).
- **GQ2 — where the editor lives.** (a) a new centre tab "Geometry"; (b) an edit mode of the Device view.
  **Recommended (a)**: the Device view shows results on the mesh, the editor shapes before any mesh; ruled (a).
- **GQ3 — 3D.** (a) a later unit extruding the 2D cross-section; (b) in this unit. **Recommended (a)**, ruled (a).
- **GQ4 — tensor decks.** (a) edited too, rectangles only (a tensor region takes boxes), the grid lines kept;
  (b) converted to `generate2d` on the first edit. **Recommended (a)** with (b) offered by a button; ruled (a).

**Slice 0 results (2026-11-20; GD1, GD2):** `mesh::layered_partition2d` (`src/mesh/layered2d.cpp`) and the deck's
`layered = on` / `domain` keys with region `box` and `polygon` layered in deck order.
- **GE2 met:** 1000 seeded random rectilinear layouts (1–30 shapes, past the domain and overlapping, five regions) tile
  exactly: simple counter-clockwise polygons, no vertex inside another polygon's edge, every unit cell's centre in
  exactly one polygon of the topmost shape's region, the areas per region exact. A ring (a hole) and a C whose ends
  touch at a corner are cut into simple pieces. A layered deck (a junction under an oxide, with a T-junction) builds
  to the same mesh digest and node regions as its hand tiling.
- **GE5 (partition part) met:** 200 random rectangles on a 1000² lattice resolve in 1.8–2.0 ms (median of 5, three
  runs; 43 polygons). The first version took 16 ms; what fixed it: each row painted top shape down with every cell
  written once, parts labelled and their boundary edges found in one pass each, span flood fills, and the breakpoint
  lines hidden under later shapes merged away first (exact: such a line carries no polygon edge).
- **GE9 for slice 0 met:** 10 hand mutants caught: the replace order reversed; T-junctions not inserted; a hole not
  cut; a pinch not cut; hidden rows judged by one column; the background ignored; the clamp past the domain wrong;
  slanted edges accepted; the deck's domain read wrongly; `domain` accepted without `layered`.
- **GE1:** the mesh suite passes (147 cases). The deck build suite's four failures are the working tree's own: three
  FinFET decks deleted from `decks/` and `finfet_3d.deck` edited (uncommitted, not this unit's), so their tests
  cannot open or match them.
- **Limitation recorded:** a cut (for a hole or a pinch) is a same-region internal edge; where one ends on the outer
  boundary it splits that boundary edge, so a patch across that point is refused by the mesher, as on a hand tiling.

**Slice 1 results (2026-11-20; GD3–GD5):** the editor's model `app::geo` (`include/NiTCAD/app/geometry.hpp`,
`src/app/geometry.cpp`; plain C++, no Win32, GE10): the deck read into regions (each a layer, its box and polygon
lines, the background the one without shapes), contacts (a generated deck's `[mesh] patch` segment, a tensor deck's
face clipped by its box) and the domain; written back by minimal edits; 17 edits; an undo history of deck texts;
the grid, snapping, the Shift constraint, handle drags and hit tests. Coordinates are integers in 0.1 nm.
- **Writing:** region sections keep their slots (layer order permutes them); each changed entry is rewritten in
  place, comments and blanks between entries kept; a removed section goes with the blank line before it; new
  regions follow the last region, new contacts end the file; the deck's line ends (LF or CRLF) are kept. A new
  contact's patch name is made unique among the deck's patches. The first edit of a tiled deck adds `layered = on`
  (refused while a polygon has a slanted edge); a polygon keeps axis-parallel edges (a dragged vertex's neighbours
  follow along their shared edges; inserted and removed vertices must lie on an edge) until slice 4.
- **GE3 met:** `mosfet_2d.deck` (tensor, CRLF), a layered deck with comments and a continued polygon, and a tiled
  deck (CRLF, last line without its end) come back byte-identical with no edit; after edits only geometry lines
  change (a line-multiset diff), and the edited layered and tiled decks build.
- **GE4 met:** 200 seeded random edits (16 kinds) on a layered and on a tensor deck: after each, the written deck
  reads back to the edited model and only geometry lines changed (a whole removed section excepted); undone one by
  one, each step gives the earlier text byte for byte, ending at the original; redo replays to the last.
- **GE5 (write part) met:** a write takes 0.03–0.04 ms mean over the 200 edits; on a 300-region deck a write and a
  read together take 0.5–0.8 ms per edit (three runs).
- **GE6 met:** at 96, 144 and 192 dpi and zooms 1e-2, 0.1, 1, 10, 100, 1e3 and 1e4 px/µm: box handles (picked within
  4.5 px of the centre, scaled; outside at 5 px), polygon vertices and edges (4 px), contact ends before the contact
  (5 px), the topmost shape, the background inside the domain, nothing outside it; model points round-trip through
  the screen. The grid step is the smallest 1, 2 or 5 × 10ⁿ µm at least 12 px apart (scaled), never below 0.1 nm.
- **Kernel change (slice 0's limitation lifted for holes):** a hole is now cut by two rays along grid lines from two
  of its corners to the part's outer loop, chosen to end on the domain's boundary off every patch's inside (a
  corner first), else on an interface, a patch's inside last; then the shorter. A straight cut across the part
  put a vertex inside a full-width body contact, and one ending on an interface made `generate2d` refuse a corner
  zone (P4b), both seen while writing the tests. New deck test: a buried box under a full-width body contact
  builds. Pinches are still cut along a column line. The mesh suite passes (122 cases, GE2's 1000 layouts included).
- **GE9 for slice 1 met:** 12 hand mutants caught (a 13th's edit did not apply and was not run): lengths not
  trimmed of trailing zeros; nm read as Å; comments between entries dropped; a removed section's blank kept; a
  vertex's neighbours not following; a tiled deck not made layered; undo not feeding redo; the snap radius not
  scaled by dpi; handles picked at 6 px; grid steps without 2 × 10ⁿ; a tensor contact's box clipped wrongly; the
  ray choice ignoring patches.
- **Found on the way:** slice 0 left `layered` and `domain` out of the app's deck vocabulary (XG2's test caught it
  once the app suite ran); added. The app suite passes (88 cases).

**Slice 2 results (2026-11-20; GD6–GD8):** the Geometry tab (`src/app/geometry_page.{hpp,cpp}`), second after Start
(Ctrl+2; the other tabs' keys move up one), on Direct2D (an `ID2D1HwndRenderTarget`, software when the window runs on
WARP) with DirectWrite labels.
- **Canvas (GD6):** rulers in µm (labels at the 1, 2, 5 step ≥ 60 px apart), the grid, the background and each
  region's shapes filled by material in layer order, outlines, the domain, contacts as thick segments with their
  names, the selection (a box's eight handles, a polygon's vertices; the background's handles are the domain's in a
  layered deck), rubber-band selection (Ctrl adds), wheel zoom at the cursor, pan by middle drag, Space+drag or the
  Pan tool, Fit (F, Home). The cursor's position in µm is in the status bar's cursor part. 1D and unknown decks show
  why; 3D decks "3D geometry editing is a later unit"; a deck without `[mesh]` offers **New 2D layout** (a layered
  template).
- **Tools (GD7):** a strip of tool buttons and Edit ▸ Geometry tools (keys V R P C H M while the canvas has the focus;
  no accelerators, so the deck editor keeps its letters). Polygon points are kept horizontal or vertical from the last
  (slice 4 brings slanted edges); Enter, a double-click or a click on the first point closes it. Double-click on a
  selected polygon's edge inserts a vertex, on a vertex removes it. Arrows nudge by a grid step (Shift ×10); Delete;
  Escape cancels or deselects. New Edit commands: Redo (Ctrl+Y), Duplicate (Ctrl+D), Delete (Del); on the Geometry tab
  Undo, Redo, Cut, Copy, Paste, Duplicate and Delete act on the canvas. Every finished gesture writes the deck at once
  (minimal edits, then the deck pane's text and the live check); undo and redo step through those texts and start
  again when the deck is edited elsewhere. A drawn contact is named `contact_<n>`, kind ohmic.
- **Inspector (GD8):** the Geometry section takes the inspector's place on the Geometry tab: a line on the selection,
  name, material (the materials the deck build knows), Raise and Lower layer, x/y/width/height (a box; the domain for
  the background), the vertex list with Apply (a polygon), a contact's kind and its ends, and the material to draw
  with. Values are in µm unless typed with a unit; Enter takes them; a refused value keeps the deck as it was. A region
  selected in the tree selects it on the canvas.
- **GE8 (slice 2's part) met:** a window test drives the canvas by messages: a rectangle drawn, moved right and down,
  resized by its handle; undo and redo through those; a polygon clicked and closed with Enter, a vertex dragged (its
  neighbours follow); a contact dragged along the bottom; duplicate and delete; the deck's other lines untouched; and
  the deck as drawn builds (the live check reports nothing). The inspector test types a height, a name, a bad width
  (refused), the domain's width and a polygon's vertex list. The MOSFET drawn and run is left to slice 3 with GE7.
- **GE5 (drag part) met:** 200 shapes and 5598 vertices (199 boxes and a 4802-vertex staircase), a box dragged 60
  frames on software Direct2D: median 2.5–2.7 ms a frame, slowest 5.3–6.5 ms (three runs; limit 16 ms).
- **GE9 for slice 2 met:** 11 hand mutants caught: the undo history not recorded; R choosing the polygon tool; a move
  without its vertical part (survived the first test, which moved only sideways; a downward drag added); a handle drag
  using the wrong handle; a drawn contact's kind; inspector numbers read as nm; the domain edited as a box; a polygon's
  last point dropped on closing; Redo undoing; a contact's kind not written; a domain change leaving a contact off its
  side. New model tests: a contact's kind, the domain widened, narrowed (a contact clipped) and refused.
- **Not done here, recorded:** the Measure tool draws its own line and readout rather than reusing the probe's ruler;
  the tree links regions only (contacts are not tree items). The app suite (89 cases) and the window suite (45 cases,
  the hidden ones not run) pass.

**Slice 3 results (2026-11-20; GD9):** live validation and the mesh preview, on `app/geometry-validation` (from main
`9b5292e`; built in its own worktree, as the main checkout's uncommitted deck deletions stop CMake's configure).
- **Validation (GD9):** the window hands the Geometry tab every live check of the deck (the build's own diagnostics,
  unchanged) and the mesh it built. `geo::locate_problem` (model, no Win32) places the build's error on the layout:
  the line's shape (or its region, the background as the domain), the contact of its section, the contact of its
  `[mesh] patch` line with every contact overlapping it, the points "(x, y)" in cm the mesher names, else the domain;
  no marker outside the geometry (a doping profile). Markers are drawn in red (outline, a wide segment, a ringed
  cross, the domain's border) and the message along the canvas's foot. On the Geometry tab the message goes to the
  console's Messages (error level), and Run is off (the toolbar, the menu and the command itself, with the reason in
  the status bar) while the text shown does not build; the other tabs keep Run as before. A deck opened and not yet
  checked is checked when the tab is shown.
- **Build change (deck layer):** a `generate2d` error of one patch (the mesher's "patch ...", "patches overlap",
  "a patch end within 3 s ...", named by index) is now reported at that `patch` line rather than the `[mesh]`
  header, so the deck pane and the canvas point at it. Only invalid decks' line numbers change (GE1: the deck suite
  passes, 123 cases).
- **Mesh preview:** a Mesh preview toggle on the tool strip overlays the last check's mesh edges while that check is
  of the text shown ("Building the mesh..." or "The mesh shows when the layout builds." otherwise).
- **GE7 met:** a table of invalid layouts made by model edits (two contacts overlapping, a contact's end within 3 s of
  a junction, a box too small for the interface resolution, an unknown material on a region and on the background,
  a contact's unknown kind, a zero-height box) and a doping error: each problem's message is `check_deck`'s, its
  line and markers as expected; four successive valid edits (a moved well, a new oxide layer, a gate contact, a layer
  raised, a rename, a wider gate) each build.
- **GE8 met (slice 3's part):** a window test draws a MOSFET on New 2D layout (silicon substrate, a 50 nm oxide box,
  source and drain contacts on the sides, a gate electrode over the oxide, a body along the bottom, each named in the
  inspector), types what the canvas does not draw (the gate's metal, the pockets' Gaussians, a five-point gate
  sweep), and Run solves it to the end. Another draws a contact over the source's patch: the canvas marks both
  contacts with the build's message, Run is off (and a direct Run starts nothing), the console has the error; undone,
  it builds and the mesh shows again.
- **GE5 still met:** the drag test's median frame 2.7–3.1 ms (slowest 5.6–5.7 ms).
- **GE9 for slice 3 met:** 13 hand mutants caught (one, a shape's line ignored, survived until a zero-height box case
  was added): patch errors not located; overlapping contacts not marked; the named points read in µm; the background
  marked as a shape; the domain marked beside the points; the check not tied to the text shown; Run not gated; the
  Run command not guarded; the problem not in Messages; a stale mesh shown; a corner drag always vertical; the
  build's messages not passed to the page.
- **Slice 2 defects found by GE8 and fixed:** a contact drag started at a domain corner always took the vertical side
  (now the side of the larger offset); a new contact after a rename reused `contact_1` (now also unique among the
  deck's patches); undo or redo of a drawn contact left the selected contact past the end (a crash or `bad_alloc`,
  seen intermittently; the selection is now pruned on every text applied); a deck opened while the Geometry tab was
  shown jumped to the View tab when its build ended (now it stays); a tab change did not refresh the commands' states.
- The app suite (91 cases) and the window suite (47 cases, the hidden ones not run) pass; full suites not run.

**Slice 4 results (2026-11-21; GD10):** slanted edges, on `app/geometry-validation` after slice 3 (same worktree).
- **Kernel (mesh layer):** a layered layout with any slanted edge goes to a general planar arrangement
  (`general_partition` in layered2d.cpp); rectilinear layouts keep the grid kernel unchanged (GE5 still 1.9-2.4 ms
  for 200 rectangles, GE2's 1000 random rectilinear layouts as before). Every shape edge and domain side is split
  where it meets another: sides and touches by orient2d, exact on the input's own points, and within a tolerance of
  1e-12 times the domain's extent counted as on the line or touching, so a vertex or edge meant to lie on a slanted
  edge but rounded off it by decimal coordinates splits it there instead of leaving a sliver. A proper crossing is
  computed once (along a horizontal or vertical segment its constant coordinate kept exactly, and any crossing within
  the tolerance of such a line put on it), then shared: within the tolerance of any input vertex or earlier crossing
  it is that point. Edges are labelled side by side (own shape by direction, others by an exact point-in-polygon of
  the midpoint), each region's loops traced and split at repeated vertices, holes cut by two bridges to one loop
  (ranked, tested in order, a near-collinear bridge refused as a sliver), and a vertex kept where three or more
  boundary edges meet or its two edges come from different lines. Shapes now need three vertices or more; a
  repeated vertex is refused.
- **Model and editor:** polygons may slant: `add_polygon` takes any simple outline of non-zero area, three vertices
  or more within 10 cm of the origin (exact 64-bit cross products of the 0.1 nm quanta); `insert_vertex` only exactly
  on the edge; `remove_vertex` while the outline stays valid; `set_vertex` gains `follow` (neighbours follow along
  shared horizontal and vertical edges; off, the vertex moves alone). A tiled deck with slanted polygons is now edited
  (made layered). On the canvas the Polygon tool takes clicks as they are, Shift keeping a point horizontal or vertical
  from the last; Ctrl+drag moves a vertex alone; the closing edge is no longer forced rectilinear.
- **Tests:** a trapezoid gate over a slanted-wall STI (shared computed crossing, then meshed by `generate_mesh`); a
  slanted ring with a hole and a crossing triangle; 1000 seeded random layouts of boxes, trapezoids and triangles (the
  last 300 in decimal cm, so no coordinate is exact in binary), each checked for simple counter-clockwise polygons, no
  vertex inside another's edge, area equal to the domain's, vertices on axis lines exact, and 200 sample points each in
  exactly one polygon of the topmost shape's region (a 6000-trial run, half decimal, also passed once); a vertex rounded
  above or below a slanted edge; the model's slanted rules; a slanted tiling made layered and built; a deck with a
  slanted well built (the deck test that expected "horizontal or vertical" now expects "repeats a vertex"); a window
  test drawing a trapezoid gate, a Shift-constrained polygon and a Ctrl-dragged vertex, the deck built.
- **Speed (recorded, no target):** 200 random trapezoids resolved in 60-120 ms (median of 3, 61 polygons); the
  slowest of the 1000 random layouts 2-5 ms. A first version took 4.5 s for the 200 (every bridge candidate tested
  before ranking); fixed by testing candidates in ranked order.
- **GE9 for slice 4:** 14 hand mutants, 13 caught: the near-touch tolerance off; proper crossings by `||`; axis
  snapping off; the two bridges not required to reach one loop; the sliver refusal off; the reversal check off; the
  exact on-edge test off in `insert_vertex`; `follow` always on; the 10 cm bound off; Shift ignored (always
  constrained); Ctrl ignored; input vertices not shared; the tolerant side off. **Survived (recorded):** the x-sweep
  prune not widened by the tolerance (no test puts a vertex within 1e-12 of a segment across an exact x gap). Two
  mutants on redundant checks (zero area given reversals) led to removing the redundant check instead.
- The mesh (152 cases), deck (9 failures as expected, unchanged), app (92) and window (48) suites pass; full ctest
  suites not run.

**GE5 on CI (2026-11-21; owner "fix them").** CI on main 9b5292e and PR #110 failed the GE5 timing checks: Debug
37-40 ms for the 200-rectangle partition (limit 5) and 28 ms per edit on the 300-region deck (limit 2); Release
2.24-2.29 ms on that deck. That test timed a write and the read after it, while GE5 bounds the write, so it now times
the write alone (the read is recorded beside it): locally 0.52 ms write, 0.23 ms read. Both GE5 checks are Release
only, as KG6 (owner, 2026-11-19); Debug records the times. The limits are unchanged. `mesh.layered2d` and
`app.geometry` pass in Release; full suites not run.

### Unit 26l — 3D solid modeller (proposed 2026-11-21; on `app/solid-modeller`)

**Owner request (2026-11-21):** "Full 3D solid modeller", picked over an extrusion editor (a 2D cross-section swept
along one axis) and over editing a 3D deck's cross-section in the 3D view. The branch starts from
`app/geometry-validation` `0e9d1ed` (Unit 26k slices 3–4, PR #110 open), since it reuses the 2D editor's model,
page and validation code.

**Where things stand.** A 3D device is typed. `generate3d` takes one convex polyhedron (a `cuboid`), split by
planar interfaces into slabs; `extrude3d` takes a conforming 2D tiling of cells, swept through stations, a region
per cell and section, contacts as patch boxes. Under both lies the polyhedral-complex generator (Unit 53: non-convex
volumes, many regions, interfaces with feature lines, patches on facets), which takes vertices, planar facets and
closed volumes, but no deck writes such a complex directly. Nothing in 3D can be drawn; the 3D view only shows
results.

**The model the modeller uses** is the 2D editor's, one dimension up. Solids are layered: a background region fills
a box domain, and every solid after it replaces what it covers. The complex the generator takes (conforming facets,
bit-identical shared vertices, edges in exactly two facets per volume) is computed from the layers, never built by
hand.

**Standing conditions:** as Unit 26k. The deck is the only source of truth: the modeller reads the deck and writes
deck text through the deck editor's undo. Every deck that builds today builds to the same bits (the new keys are
opt-in). Win32 + D3D12 + Direct2D only, no new dependency. The model has no Win32 type. Criteria are fixed before
code. Targeted Release tests per slice, full suites and Debug only on the owner's word. No screen automation.
Nothing pushed or merged without authorization.

**Left out (later units):** curved solids (cylinders, spheres; faceted, they would mostly break the generator's
60-degree interface rules); arbitrary oblique planes (see SQ4); boolean subtraction other than layering (a void is a
solid of the background region or of `void`); drawing doping; process emulation (etch, deposit); importing
STEP/GDS; relaxing the 3D generator's interface rules (mesh work, its own unit).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| SD1 | 0 | **Layered 3D decks.** `[mesh] type = generate3d` gains `layered = on` (default off). Layered, the domain is `[mesh] domain = x0 x1 y0 y1 z0 z1`. One region without solids is the background. Regions give `box = x0 x1 y0 y1 z0 z1` or `prism = <axis> <lo> <hi>, u0 v0, u1 v1, ...` (a simple polygon in the plane across the axis, extruded from lo to hi). Solids apply in deck order, each replacing what it covers. Contacts and boundary patches as `extrude3d`'s: `patch <name> = x0 y0 z0  x1 y1 z1` (boundary facets inside the box take its name). | Opt-in; meshed by the polyhedral generator. |
| SD2 | 0 | **The 3D partition kernel** (`mesh::layered_partition3d`, mesh layer, exact for axis-parallel faces): the layers resolved into a `PolyhedralGeometry` that passes `check_complex`. Rectilinear solids in slice 0 (boxes, and prisms of rectilinear polygons): an overlay of every solid's breakpoints on the three axes, each cell to its topmost solid, cells merged into maximal blocks per region (rows, then slabs). Facets split where neighbours' faces meet (T-junction vertices inserted on shared edges, so facets meet edge to edge). Interfaces and boundary facets on one plane are merged where the generator's piece rules want one polygon. | Coordinates are the inputs' own doubles; shared vertices bit-identical by construction. |
| SD3 | 1 | **The solid model** (`app::solid`, no Win32): the domain, the solids (name, material, box or prism, order), the contacts (patch boxes on boundary faces), read from a layered `generate3d` deck and written back by minimal edits (only the `[mesh]` geometry keys and the regions' `box`/`prism` lines change; comments, order and every other section kept). An `extrude3d` deck reads as one prism per cell and section run, exact (SQ5). Edits as operations: add, delete, move, resize a face, set a prism vertex, insert/remove a vertex, rename, set material, reorder, duplicate, mirror, add/move/resize a contact; each one deck edit, one undo step. Snapping: a grid on the work plane (1, 2, 5 × 10ⁿ µm from the zoom), vertices, edge midpoints, edges, faces; coordinates quantized to 0.1 nm. | Reuses 26k's quantum, deck writer and undo. |
| SD4 | 2 | **The 3D modeller** on the Geometry tab for 3D decks (SQ3). A D3D12 viewport through the render layer: solids shaded by material colour in layer order, edges outlined, the selection highlighted, a translucent mode and a section plane to see buried solids. A work plane (XY, YZ, ZX, or a picked face's plane) with its grid. Orbit, pan and zoom as the 3D view; the axes triad and view presets (26h). Picking by exact ray tests in model coordinates. A move gizmo (axis arrows, plane squares), face handles (a face dragged along its normal resizes), vertex handles on prisms. Tools: Select, Box (a rectangle on the work plane, then its height), Prism (a polygon on the work plane, then its height), Contact (a rectangle on a boundary face), Pan, Measure (3D distance). The inspector's Geometry section with numeric fields (x/y/z ranges, the prism's axis, extent and vertex list). | Hit tests in screen pixels (handles 7 px at 96 dpi, scaled). |
| SD5 | 3 | **Live validation and mesh preview:** the layout checked as it changes with the build's own checks (the partition, `check_complex`, the generator's interface and patch rules), each problem a marker in the viewport (the solid, the contact, the points the mesher names in (x, y, z)) and a line in the console's Messages; Run disabled with the reason while invalid. **Mesh preview:** the surface mesh (boundary and interface triangles) overlaid in the viewport. | As GD9. |
| SD6 | 4 | **Slanted prism faces** (trapezoidal fins, STI sidewalls): a prism's polygon may slant (26k's 2D kernel on each cross-section), so its oblique faces are parallel to the prism's axis (SQ4). | Exactness argued in SQ4. |
| SD7 | — | **Commits:** one branch, this proposal, a commit per slice, then ARCHITECTURE (the app's chapter, the mesh chapter's partition, and the 11.x row) and the records. The owner may stop after any slice. | |

**Criteria (fixed before code; a target missed is recorded as a limitation, never loosened):**

| # | Criterion |
|---|---|
| SE1 | **Nothing else moves:** every example and test deck builds to the same bits as before (the existing deck and mesh tests unchanged and passing). |
| SE2 | **Partition exact:** for 1000 seeded random rectilinear layouts (up to 20 solids, boxes and rectilinear prisms on all three axes, past the domain and overlapping, voids included), the kernel's complex passes `check_complex`. The region at each of 200 random sample points equals the topmost covering solid's (exact side tests). The volume per region equals the exact layered volume (integer lattice, so exact in double). A layered deck describing an existing `extrude3d` device (the reduced FinFET) builds, with region volumes and contact areas equal to the `extrude3d` deck's to 1e-12. |
| SE3 | **Round trip:** every 3D layered fixture deck read into the model and written back with no edit gives byte-identical text. After an edit, only the geometry lines differ (line diff), and the deck builds. |
| SE4 | **Undo:** a seeded sequence of 200 random edits, each undone in reverse, restores the deck text byte-identical at every step; redo replays it. |
| SE5 | **Speed:** an orbit or a drag (move or face resize) redraws in ≤ 16 ms per frame median on WARP with 100 solids. The partition of 100 random boxes takes ≤ 50 ms. The deck written back per edit ≤ 2 ms. |
| SE6 | **Picking:** solids, faces, edges, vertices, gizmo handles and contacts picked correctly at 96, 144 and 192 dpi, at zooms 1e-2 to 1e4 px/µm and at six camera directions (the presets) plus two oblique ones (table test). |
| SE7 | **Validation agrees with the build:** for a table of invalid layouts the modeller's message equals the deck build's error and its markers are where expected; every valid one builds. |
| SE8 | **Interaction:** window tests drive the viewport by messages (no screen automation): draw a box on the work plane and give its height, move it with the gizmo, resize a face, draw a prism and edit a vertex, draw a contact on a face, undo/redo; the deck text and device are as expected. A drawn 3D MOSFET (substrate, oxide, gate, source, drain, body; coarse, at most 20k nodes) runs its equilibrium and one gate bias. |
| SE9 | **Mutants:** at least 8 per slice caught (e.g. the replace order reversed; a T-junction vertex not inserted; a facet left unsplit; a ray test off by a face; the quantum dropped; undo of a reorder not restoring the order). |
| SE10 | **Layering:** the model in `app` without Win32 headers; the kernel in `mesh` without `app` or `deck`; no new dependency. |

**Questions (each with a recommendation; ruled by the standing "choose recommended" rule, 2026-11-11):**

- **SQ1 — what the deck stores.** (a) The layered solids (`layered = on`, `box`/`prism` lines), the user's intent
  kept, meshed through the polyhedral generator. (b) The computed complex (vertices, facets, volumes) written into
  the deck: exact but unreadable, and a reopened deck has lost which solid was drawn over which.
  **Recommended (a)**, ruled (a).
- **SQ2 — primitives.** (a) Boxes and prisms (a polygon extruded along any axis), which cover planar MOSFETs, FinFETs,
  nanosheet stacks, STI and wrap gates. (b) Also faceted cylinders and spheres. (c) General polyhedra.
  **Recommended (a)**: (b) mostly fails the generator's 60-degree interface rules, and (c) needs a general 3D
  arrangement with oblique planes (SQ4). Ruled (a).
- **SQ3 — where the modeller lives.** (a) The Geometry tab shows the 2D editor for 2D decks and the 3D modeller for
  3D decks; New 3D layout starts one. (b) A separate tab. **Recommended (a)**: one place to shape a device; ruled (a).
- **SQ4 — slanted faces.** (a) Slice 4 allows only oblique faces parallel to a prism's axis. A cut across the axis
  keeps the face's vertices exact (the cut's coordinate exact, the others the polygon's own). A cut along the axis
  meets the face in a line through a rounded 2D crossing, shared as in 26k slice 4. The split face then lies on two
  exactly planar quads meeting at a fold of about an ulp, which `check_complex` accepts. Whether the generator
  accepts such near-coplanar interface pieces is not known: where it refuses, the editor says so and the layout is
  marked (recorded, not loosened). (b) Arbitrary oblique planes: rounded points off their planes break exact
  planarity, so this needs exact rational vertices or a snap-rounding scheme, which is a research item.
  **Recommended (a)**, ruled (a).
- **SQ5 — `extrude3d` decks.** (a) Read into the model as prisms, exact (each cell's run of equal sections is one
  prism along the axis). The first edit writes the deck as layered `generate3d`, offered by a button with the change
  of mesh type said. (b) Edited in place as an extrusion. **Recommended (a)**: one 3D model; the extrusion's swept
  mesh is kept until the button is used. Ruled (a).
- **SQ6 — the generator's interface rules** (60-degree sectors, 4 s spacing, no holed pieces, no trace at a domain
  corner, …). (a) Kept: the modeller shows each refusal on the layout. (b) Relaxed in this unit.
  **Recommended (a)**: relaxing them is mesh work with its own criteria. Many drawn layouts will need a small gap or
  an overlap to pass, and SE8's MOSFET is drawn to pass. Ruled (a).

**Risks (recorded up front):** the 3D generator's rules are strict (SQ6), so whether a typical drawn device builds
depends on them; slice 3's validation is how this shows. 3D runs are slow (the FinFET takes 18 min on 6 threads),
so SE8's device is coarse. The D3D12 viewport's picking and gizmos are new code in the render path, so SE5 and SE6
are the gates there.

**Slice 0 results (2026-11-21; SD1, SD2):** `mesh::layered_partition3d` (`src/mesh/layered3d.cpp`, declared in
`complex.hpp`) and the deck's layered `generate3d`. Run under the owner's 4-core cap: builds `-j 4`, tests on 4 P
cores.
- **Kernel:**
  - **Grid:** each axis's breakpoints come from the solids and the patch boxes, merged within 1e-12 of the extent
    (20nm and 0.02um are one plane).
  - **Rows:** each row along x is resolved top solid first as run-length intervals. A box gives its index range; a
    prism along x gives its cell centre's side of the polygon; a prism along y or z gives a scanline across its
    rectilinear polygon. All of this is on grid indices, so exact.
  - **Blocks:** blocks are the identical runs of a row, extended over the next rows and slabs; each block is a volume.
  - **Facets:** the faces of the blocks on either side of each plane are overlaid (on the domain's boundary, with
    the patch boxes) and merged into rectangles. Every vertex on a facet's edge is inserted in it.
  - **Names:** boundary facets take the first patch box holding them (to 1e-12), else their face's name (x_min ...
    z_max). A patch box that selects nothing is an error.
  - **Validation:** the output passes `check_complex` by construction. The kernel no longer runs it: it took 3.2 s on
    3500 facets, the generator applies the same rules when it indexes the complex, and the tests run it on every
    output.
- **Deck (SD1):** `[mesh] type = generate3d`, `layered = on`, `domain = x0 x1 y0 y1 z0 z1` (required), and
  `patch <name> = x0 y0 z0 x1 y1 z1`. Regions give `box = x0 x1 y0 y1 z0 z1` and `prism = <x|y|z> <from> <to>, u v,
  ...`; one region without a solid is the background. Errors are placed on the solid's or patch box's line. `domain`
  without layers, `cuboid`/`vertex`/`face`/`plane` with them, and `prism` in any other mesh type are refused with a
  message. The app's deck vocabulary knows `prism`, and `layered`, `domain` and the box `patch` for `generate3d`.
- **SE2 met (kernel part):**
  - 1000 seeded random layouts (1–20 solids: boxes, and L- and U-shaped prisms on all three axes, past the domain,
    five regions, patch boxes in a third of them) give complexes that pass `check_complex`. Each volume is a box; the
    volumes fill the domain; each of 200 sample points lies in exactly one volume, of the topmost solid's region; the
    volumes per region equal the unit-cell counts exactly.
  - The reduced FinFET drawn as five boxes over a background has the same region volumes and the same source, drain,
    gate and substrate areas as `extrude`'s complex of `decks/finfet_3d_reduced.deck`, to 1e-12 (the rest of the
    boundary matches in total).
  - **Not run:** the layered FinFET deck's own build. It meshes for about 20 minutes, like the extruded one, so it
    is a hidden `[.benchmark]` case.
- **SE5 met (partition part):** 100 random boxes on a 1000³ lattice resolve in 37–44 ms (median of 3; 576 volumes,
  3435 facets). The first version took 3.3 s: 3.2 s in `check_complex`, and about 100 ms painting and merging a full
  8-million-cell grid. The run-length rows fixed the rest.
- **Deck tests:**
  - A layered MOS stack (silicon, an oxide layer, a gate patch box on the top, a body on the bottom) builds through
    the polyhedral generator, with the unnamed faces named.
  - The oxide written as a box, as a prism along z and as a prism along x meshes to the same mesh digest.
  - Eight refusals are checked with their lines.
- **SE1:** the deck suite passes (9 failures marked as expected, unchanged), the mesh suite passes (160 cases), and the
  app suite passes (92 cases, the vocabulary test included). Full ctest suites not run.
- **GE9-style mutants for slice 0 (SE9):** 14 run; 13 caught:
  - solids painted bottom up;
  - vertices on a facet's edge not inserted;
  - patch-box matching without its tolerance;
  - near-equal breakpoints not merged;
  - the prism scanline's crossings inverted;
  - a prism's (u, v) swapped;
  - blocks merged across regions;
  - the faces' default names swapped;
  - an unused patch box allowed;
  - the slanted-edge refusal dropped;
  - errors placed at the section's line;
  - `layered = on` ignored;
  - every prism along z.

  Two of them (the patch tolerance, the breakpoint merge) survived until a test wrote one plane two ways. One is
  equivalent and recorded: the solids not sorted by line (a region's solids are contiguous in the deck, so they
  already come in line order).

**Slice 1 results (2026-11-21; SD3):** `app::solid` (`include/NiTCAD/app/solid.hpp`, `src/app/solid.cpp`), the 2D
editor's model one dimension up. Run under the owner's 4-core cap: builds `-j 4`, tests on 4 P cores.
- **Model:** the domain, regions with their solids (a box, or a prism: axis, extent, polygon in (u, v)), and contacts
  as patch boxes. Coordinates are 26k's 0.1 nm quanta, and undo is 26k's `EditHistory`. The 2D writer's line edits
  moved into a shared internal header (`src/app/deck_text.hpp`), so both models now use them. The 2D model's behaviour
  is unchanged (its tests pass). `geo::polygon_ok` is now public for the prisms.
- **Reading and writing:** a layered `generate3d` deck is written back by minimal edits:
  - **What is rewritten:** region `box`/`prism` lines, headers and materials; the `[mesh]` `layered`, `domain` and
    patch lines; the contacts' sections.
  - **Patch boxes:** a patch box no contact names is kept. One whose contact was deleted is dropped.
  - **Decks not drawn here:** a tensor 3D deck, an unlayered `generate3d` deck and 2D decks read as unsupported, each
    with its reason.
- **`extrude3d` (SQ5):**
  - **Reading:** read exactly as one prism per cell and `extrude` line, and not edited in place.
  - **Conversion:** `convert_to_layered` gives the note for the button. Writing the converted model rewrites the
    following:
    - `type` becomes `generate3d`, with `layered = on` and the domain added;
    - `axis`, `stations`, `cell` and `sweep` are dropped, and a swept deck's `lines`, `straddle`, `layer`,
      `min_angle` and `box` with them;
    - each `extrude` line becomes prisms.
  - **Refusals:** conversion is refused with the reason when a cell is void in a section, when a region holds no
    cell, or when the cells do not fill their bounds.
  - **Kernel check:** the reduced FinFET converted (with its C-shaped cells and collinear vertices) passes the
    slice 0 kernel and `check_complex`, with its volume per region exact.
- **Edits:** each is one deck write and one undo step:
  - add a box or prism;
  - move;
  - move a face: a box's six faces, or a prism's caps and its sides along their normals;
  - set, insert or remove a prism vertex (neighbours follow);
  - delete, rename, set the material, raise and lower;
  - duplicate, or mirror across an axis;
  - add, set, move, delete, rename a contact or set its kind (a patch box that reaches a domain face with area);
  - set the domain (contacts on a face move with it).

  A result that is empty, outside ±10 cm or slanted is refused (slanted faces are slice 4). An extruded model refuses
  edits until it is converted.
- **Snapping:** the order is a vertex (the solids', the contacts' boxes', the domain's), then an edge's midpoint, then
  the nearest point of an edge, then a face (a box or domain face, a prism cap or axis-parallel side). Faces in the work
  plane's own plane are left out. After faces comes the work plane's grid (26k's step), then the point as given. The
  tolerance is in quanta; the viewport converts 8 px at the point's depth (slice 2).
- **SE3 met:**
  - Eight decks are written back byte-identical with no edit: the well and MOS fixtures, an extruded stack, a CRLF copy,
    and the four FinFET extrusions.
  - After edits, only geometry lines differ, and a read of the written text gives the model again.
  - The edited MOS stack builds, and the converted coarse extrusion builds.
- **SE4 met:** 200 seeded random edits on the well deck, and 200 on the converted reduced FinFET (the conversion its
  first step), are each undone byte-identical at every step and redone to the end.
- **SE5 met (write part):**
  - The mean write per edit is 0.05–0.07 ms (well deck) and 0.10–0.15 ms (FinFET).
  - A deck of 100 solids is written and read again in 0.34 ms per edit.
  - The frame and partition parts are slice 2's and slice 0's.
- **SE10 met:** the model's files name no Win32 or D3D12 header or type (a test reads them), and `app_layering`
  passes.
- **Suites:** the app test executable passes 101 cases (the 2D model's included). The deck and mesh suites are
  untouched. Full ctest suites were not run.
- **Mutants (SE9):** 14 run, 14 caught:
  - a changed solid not rewritten;
  - the box text's axes out of order;
  - (u, v) swapped;
  - the void check inverted;
  - an unused patch box dropped;
  - neighbours not following a vertex;
  - a prism side moved along the wrong coordinate;
  - the work plane's faces kept in snapping;
  - a contact with area on one axis only;
  - the domain's low faces not carrying their contacts;
  - cells kept on conversion;
  - raise without the swap;
  - an extrude range reversed;
  - slanted edges accepted.

  One mutant survived until a test moved the domain's low face: the domain's low faces not carrying their contacts.
  One first failed to build and was rerun as a compiling equivalent: neighbours not following a vertex.

**Slice 2 results (2026-11-21; SD4):** the 3D modeller on the Geometry tab. Run under the owner's 4-core cap.
- **Files:** the view without Win32 is `include/NiTCAD/app/solid_view.hpp` and `src/app/solid_view.cpp`: the scene, the
  handles, picking and the drags' geometry. The window is `src/app/solid_page.hpp` and `src/app/solid_page.cpp`.
  `render` gains `view_frame`, `screen_ray` and `um_per_pixel`, the projection `draw()` uses, in double precision.
- **Which page:** the Geometry tab shows the modeller for a layered `generate3d` or an `extrude3d` deck, and the 2D
  editor otherwise. A tensor 3D deck keeps the 2D page's message. The Edit menu's undo, redo, duplicate, delete and
  tool commands go to whichever page is shown. Run gating is still the 2D page's (live validation is slice 3).
- **Drawing:** every solid is drawn whole, shaded by material. Each layer's faces are pushed out by 1e-5 of the domain
  per layer, so the later layer shows where faces meet, as the deck paints them. Edges are outlined. Contacts are drawn
  on the boundary, the selection in orange on top. Also drawn: the work plane's grid, a see-through mode (T; the
  selection stays opaque), a section plane (the cut faces are left open) and the axis triad. The camera is the 3D
  view's: right drag orbits, middle drag pans, the wheel zooms at the cursor, F fits.
- **Tools:** Select (gizmo arrows and plane squares, face handles, prism vertex handles; Shift adds, Ctrl drags a
  vertex alone), Box, Prism (edges kept along u or v), Contact (on a domain face), Pan and Measure. Arrows and
  Page Up/Down nudge. The inspector has the ranges, the prism's axis, extent and vertices, the material, raise/lower,
  the contact's kind, the drawing material, the work plane, the section and Convert. An `extrude3d` deck is shown and
  can be selected; an edit is refused with the reason until Convert.
- **SE5 met (frame part):** with 100 solids on WARP, measured over 31 frames each, the scene build included:
  - run alone (five runs): orbit median 10.2–10.8 ms, drag median 10.3–12.1 ms;
  - inside the whole window test executable (four runs): 14.1–14.4 ms and 13.5–14.5 ms. One further run read 16.7 and
    17.3 ms. That was a loaded run: the same run failed the existing hover timing check, which passed on the rerun.
  The scene build is 0.4–0.55 ms of that. The rest is the WARP present.
- **SE6 met:** a table over 3 dpi × 4 zooms (1e-2, 1, 1e2 and 1e4 px/µm) × 8 views (the six presets, the isometric and
  one more oblique view) picks, on the face toward the viewer:
  - a box's face, a prism's cap, the notch of the prism's L (the background: the exact polygon test), a contact and
    the background, and nothing off the domain;
  - past a section plane, the far faces (axis-aligned views);
  - every vertex, face, arrow and plane handle, and the edges, where no earlier-ranked handle lies near on the screen.
  Counts are checked: at least 96 each of faces, contacts, vertices, arrows and edges. The picking frame matches
  `render::screen_position` within 0.02 px in 8 views, orthographic and perspective. `on_plane` and `along_axis` give
  the point back within 1e-7 µm. Triangulation is exact on L, U, staircase and collinear polygons.
- **SE8 met (drawing part):** window tests drive the view by messages, with no screen automation:
  - a box drawn on a work plane and given its height, moved 0.2 µm by the gizmo arrow, resized by a face handle;
  - an L prism clicked and closed by Enter, with its height; a vertex dragged, its neighbour following;
  - a contact dragged on the y_max face;
  - all six edits undone to the opened text byte for byte, and redone from the Edit menu.
  The deck text after each edit is as expected, and the final deck builds to a 3D device with 4 regions and 3
  contacts. Further tests cover the inspector's fields, the material, a refused value, the work plane, Delete, and an
  `extrude3d` deck (an edit refused with the reason, then Convert and undo). The drawn 3D MOSFET run is slice 3's.
- **Found while testing:** the mesher (`complex_interfaces.cpp`) has two limits a drawn layout can hit:
  - **Interface holes:** it refuses an interface piece with a hole. A box set inside a layer's interface plane, not
    reaching the boundary, punches one into that plane.
  - **Short traces:** it refuses a trace shorter than 6 × `interface_resolution`.

  The modeller lets both through. Slice 3's live validation will mark them with the build's own checks (SD5). The
  test's box therefore starts at the oxide's corner, and its fixture uses 10 nm interface resolution.

  The edge snap keeps the raw position along an edge, as the 2D editor's does, so a click on an edge can land off the
  grid (0.1002 µm). Left as is, for consistency with 2D.
- **SE10 met:** the view files name no Win32 or D3D12 type (the slice 1 check now reads them too), and `app_layering`
  passes.
- **Suites:** the app test executable passes 104 cases, and the app window test executable 52. Full ctest suites were
  not run.
- **Mutants (SE9):** 15 run, 14 caught:
  - the layer tie rule reversed;
  - the point-in-polygon parity forced;
  - the handle radius widened;
  - the screen's y flipped;
  - `along_axis`'s sign;
  - `on_plane`'s sign;
  - the section's side;
  - the edge radius zeroed;
  - the height ignored;
  - the vertex follow inverted;
  - the face drag's sign;
  - the page chosen for layered decks only;
  - undo history not recorded;
  - the extrusion's edit guard removed.

  The survivor is equivalent: it flipped the crossing test's direction, and both directions of a closed polygon give
  the same parity. It was replaced by the forced-parity mutant. The edit guard's mutant survived at first: the writer
  also leaves an unconverted extrusion unchanged. A check on the refusal message now catches it.

**Slice 3 results (2026-11-21; SD5):** live validation and the mesh preview. Run first on the 4-core cap; the owner
lifted it to 6 cores during the slice, and the mutants and final runs used 6.
- **Where an interface fails:** the generator's interface and patch-edge rules name a facet of the complex (its
  index). For a layered `generate3d` deck the build now ends such a message with `near (x, y, z)`: the facet's vertex
  centroid in cm. Only that path changes; other decks' messages and every mesh are as before.
- **Problems placed (`solid::locate_problem`, no Win32):** the 2D editor's rules, one dimension up:
  - a region's `box` or `prism` line marks that solid, another line of the region all its solids, and the
    background the domain;
  - a contact's section marks its contact, and a `[mesh]` patch line the contacts using it;
  - a named point is marked, with the solids whose closed bounds hold it;
  - another `[mesh]` line marks the domain. Nothing is marked outside the geometry.
- **On the page:** the problems' markers are drawn in red over everything. A patch box on no boundary face is
  outlined as a box. The first problem heads the inspector with "Run is off until the layout builds"; Run is
  disabled with that reason, and the console's Messages get the line (the main window's existing path). The page
  shows the check state of the text shown (pending, valid, invalid), as the 2D page does. The tab starts a check when
  it opens on an unchecked deck.
- **Mesh preview:** a Mesh button on the tool strip. The live check's generated tets give the surface: the faces on
  the boundary and those between tets of different regions. Each tet's region is where the model puts its centroid
  (`solid::region_at`). The preview shows only for the text it was built from. Interface nodes come in pairs either
  side of the interface (FD7), so the interface faces straddle it rather than lie in it.
- **SE7 met:** a window-test table through the live check covers:
  - two valid layouts (the stack, and a box reaching the boundary), each checked valid and built with Run on;
  - four invalid ones: a flat box (its solid marked), a patch box on no boundary face (its contact), a degenerate
    domain (the domain), and a box inside the oxide's interface plane (a hole: the named point on z = 0.4 µm, and the
    oxide whose face it is).

  For each, the page's message equals `check_deck`'s byte for byte, the markers are as listed and drawn, Run is off
  with the reason, and asking for Run starts nothing. A model test covers `locate_problem` (each marker kind) and
  `region_at` (the later solid wins an overlap; boundary points count as outside).
- **SE8 met (run part):** a 3D MOSFET drawn on the page runs. The substrate is the opened deck's background. The
  oxide is drawn as a box on a work plane, given its height, and has its bottom typed as 0.55 µm. Source, drain, gate
  (set to electrode) and body are dragged on the domain's faces, each seen straight on, and named. The gate's metal,
  the source and drain pockets' doping and the sweep are typed, as 26k's GE8 does. The check is valid and Run is on.
  The mesh has 3930 nodes (limit 20k). It runs the equilibrium and the gate at 0 and 0.5 V, all converged, in
  5.8 s on 6 threads.
- **Mesh preview test:** the preview shows the built mesh, including interface edges near z = 0.4 µm away from the
  boundary. After an edit it shows "Building the mesh..." until the next check. A mesh tagged with other text is not
  shown.
- **SE1:** the deck tests for layered decks pass. The hidden `[.benchmark]` reduced-FinFET test fails before
  building: the worktree's CRLF deck defeats its text replacement. That test was not run in slice 0 either.
- **Suites:** the app test executable passes 105 cases; the app window test executable passes 55. On 6 cores SE5's
  frame medians read 5.0 ms (orbit) and 4.9 ms (drag). Full ctest suites were not run.
- **Mutants (SE9):** 12 run, all caught:
  - the solid line's index shifted;
  - the contact's section origin shifted;
  - the named point's units;
  - `region_at`'s order reversed;
  - the `[mesh]` fallback to the domain removed;
  - the check state ignoring the text;
  - Run not gated by the 3D page;
  - interface faces left out of the preview;
  - markers not drawn;
  - the build's facet point dropped;
  - the mesh tagged with the wrong text;
  - the patch-box outline removed.

  Three survived at first, until tests gained a solid overlap, an interior interface edge and a mismatched mesh text:
  `region_at`'s order, interface faces left out, and the wrong mesh tag. One first failed to build and was rerun as
  a compiling equivalent: the facet point dropped.

**Slice 4 results (2026-11-21; SD6):** slanted prism faces, on `app/solid-modeller` after slice 3 (same worktree). Run
under the owner's 6-core cap: builds `-j 6`.
- **Kernel (mesh layer):** a layered layout with a slanted prism goes to a new path in `layered3d.cpp`; layouts without
  one keep the rectilinear kernel unchanged (SE5's 100 boxes 26.7 ms, SE2's 1000 random layouts as before).
  - **One axis:** every slanted prism must share one axis A (SQ4). Slanted prisms along two axes are refused on the
    second one's line: "slanted prisms along two axes (every slanted face must be parallel to one axis, the one all
    slanted prisms share)". Boxes and rectilinear prisms on any axis stay allowed.
  - **Slabs:** the breakpoints along A cut the layout into slabs. In each slab every solid is a 2D shape in (u, v): a
    box's rectangle, a prism along A's polygon, and a rectilinear prism across A as the rectangles a scanline at the
    slab's centre gives (exact, on breakpoint indices).
  - **One arrangement for all slabs:** 26k's slanted 2D arrangement (`general_partition`) now takes a caller's label
    and pinned points (`detail::labelled_partition2d`, `src/mesh/arrangement2d.hpp`); the 2D path is unchanged. A
    point's label is its signature: its region in every slab, and the patch box holding it on each end plane. So the
    arrangement's polygons are the pieces across which nothing changes along A, simple, conforming, and each
    crossing is computed once for every slab.
  - **The complex:** Unit 53's sectioned extrusion (`extrude`) of those polygons through the breakpoints. Lateral
    facets are polygon edges swept through one slab: vertical quads, exactly planar even through a rounded crossing.
    Boundary facets take the first patch box holding them, else their face's name (x_min ... z_max), as the
    rectilinear kernel. A patch box on an end plane splits the arrangement along its edges; on a side face, its edges
    are pins on the domain's sides. `extrude` runs `check_complex`, so this path costs more (below).
  - **Coordinates:** on u and v, values within 1e-12 of the extent are first made one (20nm and 0.02um), the domain's
    ends kept exactly.
- **What SQ4 left open, answered:** a sidewall crossed by a plane (the gate's bottom across a trapezoidal fin's
  sidewall) is split at a rounded crossing. The two pieces are parallel to rounding but not exactly coplanar, and
  the generator's separation rule (kept, SQ6) measures them as two interfaces at distance zero: "two interfaces
  within 4 s". Recorded, not loosened. A sidewall no plane crosses meshes: the same fin with the gate down to its
  base gives 6517 nodes. Live validation marks the refusal at its facet as any other generator error.
- **Model and editor:**
  - `prism_ok` takes any simple polygon of non-zero area (three vertices or more).
  - `insert_vertex` takes a point exactly on any edge (64-bit cross product).
  - `move_face` still refuses a slanted side (its normal is off the quanta's lattice); its vertices move instead.
  - The Prism tool takes clicks as they are; Shift keeps the edge along u or v from the last point. Ctrl still drags
    a vertex alone. Snapping onto a slanted side face is not added (recorded).
  - **Convert:** an `extrude3d` deck with slanted cells now builds after Convert. The extrusion named only the cells
    inside a contact's box; a layered deck names the part of every face inside it. So Convert's note now names each
    contact whose box cuts across a cell's end on an end plane. The shaped FinFET's source and drain cut across the
    spacer's slanted cells, so those contacts now also take spacer slivers.
- **Tests:**
  - **Hand cases:** a trapezoidal fin under a wrapping gate (region volumes exact to 1e-12 against the hand values,
    patch areas, slanted facets present); an end-plane patch box across the fin's sidewalls and a side strip, their
    areas exact; a trapezoid, a triangle, a box corner rounded onto a slanted edge and rectilinear prisms across the
    axis (1000 sample points); one coordinate written two ways; refusals (two axes, uncovered, unused patch box).
  - **1000 seeded random layouts** (1–12 solids: boxes, rectilinear L prisms on all axes, trapezoids and triangles
    along one axis per layout; the last 300 in tenths; end and side patch boxes in a third): each passes
    `check_complex`, its volumes fill the domain, and 200 sample points each lie in exactly one volume, of the
    topmost solid's region.
  - **Generator:** the uncrossed fin meshes; the crossed one's refusal is pinned by a test.
  - **Model:** the slanted edits; the shaped FinFET converted, its note, and its region volumes through the slanted
    kernel equal to its prisms' (the extrusion's) to 1e-12.
  - **Deck:** a trapezoidal oxide (prism along y) builds. The deck test that expected "a prism's polygon has a slanted
    edge" now expects the build.
  - **Window:** a trapezoid drawn on the work plane y = 0 with a slanted click and a Shift-kept one; the polygon is
    as clicked, and the deck builds a 3D device with 3 regions.
- **Speed (recorded, no target):** the slowest of the 1000 random slanted layouts took 869 ms, most of it in
  `check_complex` inside `extrude`. The rectilinear path keeps SE5.
- **SE1:** nothing else moves. The rectilinear kernel and the 2D path are unchanged, and the layered deck tests pass.
  The hidden `[.benchmark]` FinFET case still fails before building on the CRLF worktree deck, as in slice 3.
- **Suites:** the mesh cases tagged `[layered2d]`, `[layered3d]` and `[complex]` (21) pass, the app test executable
  (105 cases) passes, and the window cases tagged `[solid]` and `[geometry]` (15) pass. Full ctest suites were not
  run, and a run of whole test executables was declined by the session's permission check.
- **SE9:** 17 mutants, all caught: the top solid first-wins; end-plane patch rects dropped; pins not made; pins not
  kept; the coordinate snap off; the scanline crossing inverted; the cross-axis rect's (u, v) swapped; a box one slab
  too many; the own-edge-by-direction rule off on the labelled path; face names swapped; the two-axes refusal off; the
  background taken for an uncovered point; `insert_vertex` off the edge; `prism_ok` rectilinear again; Convert's note
  on every meeting cell; Shift ignored; always constrained. Two first failed to build (an unused variable) and two
  ran no test (Catch2 refused the window case's full name as a filter); all four were rerun as compiling
  equivalents, the window case by a name prefix.

**Test speed (2026-11-21; owner "optimize").** The full Debug suite was first run on this branch at the owner's word
and took hours: CI's last green Debug job ran 7695 s, with `solve.quantum_confinement` 6007 s,
`adapt.adapt_derived` 2720 s, `render.image` 2577 s, `solve.energy_balance` 2430 s and `linalg` 1175 s; in Release
`solve_mosfet` 2334 s and `deck_lombardi` 1001 s. Commit `d974da3`:
- `nitcad_linalg` (Eigen's sparse LU and ILUT, AMG) is built `/O2 /Ob2` in Debug too (target property
  `NITCAD_FAST_DEBUG`); every other target keeps `/Od /Ob0 /RTC1`, and the Debug CRT, checked iterators and asserts
  stay on, so the ABI is one. Locally, Debug: quantum confinement 5928 s → 39 s, derived goals 2654 s → 24 s,
  `render.image` 48 s, energy balance 28 s, `linalg` 29 s; the whole Debug suite 835 s.
- `solve_mosfet` and `deck_lombardi` are one ctest entry per case (`catch_discover_tests`, `COST 1000`), so `ctest -j`
  spreads them; the slowest case is 392 s. CI's 4 jobs are unchanged.
- Debug only: the two 1000-layout sweeps run every 20th seed (Release all); SE5's timing limits (and the 2D
  partition's, as main's GE5 fix) are checked in Release only, as KG6 (owner, 2026-11-19), the times recorded. The
  limits are unchanged.
- A C4702 in `solid_window_test.cpp` (slice 2) broke the Debug build; fixed.
- Suites: Release 277/277 (1014 s; the PARDISO case skipped locally, no MKL runtime configured in this worktree).
  Debug: after the timing guards, all pass; `app.solid` failed once with no assertion output and passed four reruns
  (recorded, not reproduced). Merged main (PR 110) as `95b243b`; pushed; PR #111 opened (owner "push", "open pr").

**SD7 (2026-11-21):** ARCHITECTURE.md updated: 6.11 (layered 3D decks, the partition kernel, slanted prisms, the
limits found), 6.14 (as built, the modeller) and the two 11.3 rows (geometry and meshing; visualization and GUI).
Left for later units: curved solids, oblique planes, slanted sidewalls crossed by a plane (the generator's 4 s rule),
interface pieces with holes, snapping onto slanted faces.

### Unit 26m — minimal modeller tools (proposed 2026-11-21; on `app/modeller-tools`)

**Owner request (2026-11-21):** after "too many bugs, visually bad" on the 26l modeller and the drag-to-start fix,
the owner asked for a minimal, easy modeller: tool groups (Solid, Edges, Faces, Transform, and later Advanced and
Semiconductor), then "check Sentaurus, Silvaco and other CAD modellers … improve the plan", then "ignore the
emulator things" (no process emulation), then a button that fixes the camera, with the view moved by the mouse when
it is off. The design canvas the owner approved ("Its good plan it") is the private artifact
`DkqARw6W5TT8wz9nKHsGXr` (two artboards: the Geometry tab and the 23 tools). The branch starts from
`app/solid-modeller` `a8ca49f` (Unit 26l, PR #111 open), since every tool edits 26l's model.

**Where things stand.** The 3D modeller (26l) has Select, Box, Prism, Contact, Pan and Measure in a fixed tool row,
a move gizmo, face handles, prism vertex handles, duplicate and mirror as model operations, slanted prisms, live
validation and the drag-to-start gesture. The camera orbits on the right button and pans on the middle one; nothing
fixes it. The 2D editor (26k) has its own tool row. The deck already has `[doping]` profiles (`constant`,
`gaussian`, `erfc`, with axis, centre, width and lateral spread) and the mesh keys `max_size` and `layer`; it has no
named values. Solids only replace what they cover.

**What the canvas fixed, and where it came from:**
- **Selection decides the tools.** Picking a solid, an edge or a face opens its group as a small floating bar beside
  it; with nothing picked, Solid. Six groups, at most five tools showing.
- **Every tool is one drag.** The value shows at the cursor; a typed number makes it exact; Enter applies, Esc
  cancels; a refusal (the build's own check) shows before the button is released.
- **From Sentaurus Structure Editor:** an overlap rule per solid ("replaces" or "fills only empty space") in place
  of Union, Subtract and Intersect; the deck line each action writes is shown as it is dragged (as SDE echoes Scheme);
  dragged values can be named parameters (as Workbench substitutes them); doping and refinement are drawn like
  geometry.
- **Camera:** a "Camera fixed / Camera free" button; free, the mouse moves the view of the 3D model or the 2D
  layout; fixed, drags only edit the model.

**Standing conditions:** as Unit 26l. The deck is the only source of truth, written by minimal edits through the
deck editor's undo. Every deck that builds today builds to the same bits (new keys opt-in). Win32 + D3D12 + Direct2D
only, no new dependency. The models have no Win32 type. Criteria fixed before code. Targeted Release tests per
slice; full suites and Debug only on the owner's word; builds `-j 6`, runs on 6 threads. No screen automation.
Nothing pushed or merged without authorization.

**Left out:** process emulation (deposit, etch, fill, planarize, spacer, masks; owner, 2026-11-21); sweeping
geometry parameters across runs (a batch-run unit); variable fillets and rounded 3D corners (boxes and prisms cannot
hold them); oblique planes and fillets on edges not parallel to a prism's axis (26l SQ4); relaxing the 3D
generator's rules; GDS/STEP import.

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| MT1 | 0 | **Camera lock and mouse navigation** in the 3D modeller and the 2D editor. A toggle button "Camera free / Camera fixed" (lock icon, pressed when fixed; shortcut L). Free: left-drag on empty background orbits in 3D and pans in 2D; Shift+left-drag and middle-drag pan; right-drag orbits (as now); the wheel zooms about the cursor. Fixed: no mouse drag or wheel moves the view; every drag goes to the model; view presets and Reset view still work. A Reset view button. The state is per view and kept in the view settings. | Orbit pitch kept off the poles as now. |
| MT2 | 0 | **The context toolbar.** The fixed tool row becomes a group switcher (Solid, Edges, Faces, Transform, Doping, Mesh). The armed group's tools float in a small bar beside the selection's screen bounds (clamped inside the viewport, not over the picked thing); with nothing picked, in the viewport's top-left corner. Picking a solid, an edge or a face opens its group; the bar's line names the armed tool and how to drag it. Section and Measure move to the view's corner. The same switcher on the 2D editor. | Pure presentation over 26l/26k tools. |
| MT3 | 0 | **The live deck line and the value chip.** The status bar shows the deck line the current drag will write (the same text the minimal edit produces), updated each frame. The chip at the cursor shows the value; digits typed during a drag set it exactly (units as the deck's), Enter applies, Esc cancels. | |
| MT4 | 1 | **The overlap rule.** A region's solid may say `overlap = fill` (default `replace`): it takes only the cells no earlier solid holds (the background's), in both layered kernels (2D 26k, 3D 26l). Set from the inspector's two-way switch or the Solid group's Overlap tool (click a solid). | Exact: decided per cell on grid indices, as the replace rule. |
| MT5 | 1 | **Transform group:** Move (gizmo, as now), Copy (Alt+drag, or the tool: drag and release to place; one new region), Mirror (click a plane through the solid's centre or a picked face; Shift: a quarter turn about the view axis instead), Array (drag the pitch, the wheel sets the count; writes one region per copy), Align (click a face; the selected solids' matching faces move to its plane). | Each one deck edit and one undo step. |
| MT6 | 2 | **Faces group:** Push/Pull (the face handle, as now), Inset (drag inward for the margin: a new solid on the face, inset, of the chosen material, then pushed or pulled), Extrude (drag off a face: a new solid on it), Contact (as now). | Inset and Extrude are new solids, so they obey the overlap rule. |
| MT7 | 2 | **Edges group and Cylinder,** on edges parallel to a prism's axis (a box is first written as its rectilinear prism): Chamfer (one slanted segment), Fillet (an arc of segments, MQ4), Round vertex (a prism polygon's corner, or a 2D corner in the 2D editor), Taper (drag the edge sideways; the wall angle shows; 26l slanted prisms). Cylinder in the Solid group: a prism whose polygon is a regular polygon (MQ4), drawn by dragging the radius. | Exact as 26l SD6; the generator's refusals shown, not loosened. |
| MT8 | 3 | **Named parameters.** `[parameters]` holds `name = <length or number>`; a length value anywhere in `[mesh]`, region solids and `[doping]` may be `$name` or `$a + $b`, `$a - 5nm` (MQ1). Pressing N during a drag names the value (a new parameter, or an existing one picked); the model then writes the reference, and dragging a tied dimension edits the parameter's line. The inspector lists the parameters, editable. | Substituted before the build parses values, errors placed on the referring line. |
| MT9 | 4 | **Doping group:** Region profile (click a solid: a `constant` profile over its box; drag up or down for the level, a log scale), Implant (drag a window on a face, then its peak depth and lateral spread: a `gaussian` profile), Junction (drag where the junction sits: an `erfc` profile's depth). Profiles show as translucent boxes coloured by concentration; each writes one `[doping]` profile by minimal edit. | Uses the existing profile keys only. |
| MT10 | 4 | **Mesh group:** Max size (drag anywhere: `max_size`; the node count estimate shows), Interface layers (click an interface: `layer`, drag for the count and first thickness), Refine box (drag a box, then inside it for the size: a new opt-in `[mesh]` key `refine = <lo point>  <hi point>  <size>`, handed to the generators' existing `SizeBox` list, 2D and 3D; MQ5). | No generator change. |
| MT12 | 2 | **Meshable by construction** (the 3D generator's rules, in units of `interface_resolution` s: feature lines ≥ 6 s, unrelated interfaces ≥ 4 s apart, sectors ≥ 60°, no holed interface piece, no interface ending on a domain edge or corner, pieces wide enough for the 1.5 s inset). (1) **Automatic s:** the modeller computes the largest s the layout allows (the shortest feature over 6, the thinnest gap over 4, the narrowest piece over its inset) and shows it; written as `interface_resolution` unless the deck sets one. (2) **Adjusted before release:** a fillet's or cylinder's segment count drops until each segment is ≥ 6 s; Array, Align, Copy and Inset snap to keep a 4 s margin from unrelated interfaces and from the domain's edges; what cannot be adjusted (a holed piece, a sector under 60°) is shown as the refusal during the drag, naming the rule. (3) **Node estimate:** the expected node count (from s, `max_size` and the feature lines' graded pairs) shows live in the chip and the status bar. | The rules themselves unchanged here; lifting them is Unit 64. |
| MT11 | — | **Commits:** this proposal, a commit per slice, then ARCHITECTURE (6.14, the deck chapter for `overlap` and `[parameters]`, 11.3) and the records. The owner may stop after any slice. | |

**Criteria (fixed before code; a target missed is recorded as a limitation, never loosened):**

| # | Criterion |
|---|---|
| MC1 | **Nothing else moves:** every example and test deck builds to the same bits; the 26k and 26l test suites pass unchanged except where a test names the removed fixed tool row (recorded). |
| MC2 | **Camera lock:** window tests by messages: with the camera fixed, left, middle and right drags and the wheel leave the camera bit-identical and a left drag on a face still resizes it; free, a left drag on empty background changes yaw/pitch in 3D and the pan in 2D, the wheel the zoom; Reset view restores the start camera exactly; the state survives a reload of the deck. |
| MC3 | **Context toolbar:** for each selection kind the right group opens; the bar's rectangle lies inside the viewport and does not overlap the selection's screen bounds at 96, 144 and 192 dpi and six camera directions (table test, except where the selection fills the view: then the corner). |
| MC4 | **Deck line:** for every tool, the line shown during the drag equals, byte for byte, the line the applied edit writes. |
| MC5 | **Overlap rule exact:** 1000 seeded random layouts in 2D and in 3D with `fill` solids mixed in pass `check_complex` / the 2D tiling check; the region at 200 random points equals a brute-force layered resolution with the fill rule; `fill` on a solid that only covers background gives the same bits as `replace`. |
| MC6 | **Transform, Faces, Edges:** each tool's edit round-trips (undo restores the deck byte-identical; redo replays), and a table of cases builds (Array of 3 fins, Copy, Mirror, Align, Inset, Extrude, Chamfer and Fillet on a fin's top edges, Taper, a Cylinder nanowire), or is refused with the build's own message where the generator's rules refuse it (recorded). |
| MC7 | **Parameters:** decks with references build to the same bits as the same decks with the values written out; an undefined name, a cycle, or a non-length where a length is needed is an error on the referring line; changing `Lg` in the inspector moves every tied face. |
| MC8 | **Doping and mesh tools:** each writes exactly the profile or mesh line expected (table), and the built device's net doping at sample points equals the profile's formula to 1e-12 relative. |
| MC9 | **Speed:** with 100 solids on WARP, an orbit frame ≤ 16 ms median (26l SE5 unchanged); the floating bar and the deck line add ≤ 1 ms per frame. |
| MC10 | **Interaction:** a FinFET drawn by window messages only (box, Array of fins, Fillet on the fin tops, Extrude a gate over them with `fill`, Contacts, a Gaussian source/drain implant, one named `Lg`) builds and runs its equilibrium (coarse, ≤ 20k nodes). |
| MC11 | **Mutants:** at least 8 per slice caught (e.g. the fixed camera still zooming on the wheel; `fill` taking a non-background cell; a parameter substituted on the wrong line; the deck line built from stale drag state). |
| MC12 | **Layering:** models in `app` without Win32 headers; kernels in `mesh` without `app` or `deck`; parameters in `deck`; no new dependency. |
| MC13 | **Meshable by construction (MT12):** for 500 seeded random layouts drawn with the tools, every layout the modeller accepts builds, and every one it refuses is refused by the build with the same message (no false accepts; false refusals counted and recorded). The automatic s is the largest that passes: s × 1.01 is refused by the build for each layout in a table. The node estimate is within ±30 % of the built mesh's node count on that table (recorded, not a gate, if missed). `refine` boxes build to the same bits as the same `SizeBox` list given to the generator directly. |

**Questions (each with a recommendation; ruled by the standing "choose recommended" rule, 2026-11-11):**

- **MQ1 — parameter syntax.** (a) `[parameters]` with `$name` references and only `+`/`-` of lengths, substituted
  before parsing; the model writes `$x0 + $Lg` for a tied far face. (b) A full expression language (×, ÷,
  functions). **Recommended (a)**: enough for tied dimensions and readable; (b) later if sweeps need it. Ruled (a).
- **MQ2 — camera default.** (a) Free on open, fixed by the button or L. (b) Fixed on open. **Recommended (a)**: the
  owner's complaint was not being able to move the view. Ruled (a).
- **MQ3 — left-drag on empty background when free.** (a) Orbit in 3D, pan in 2D; drawing a new solid needs a tool
  armed (or an empty deck's drag-to-start, unchanged). (b) Left-drag always rubber-band selects. **Recommended (a)**,
  as the canvas; rubber-band selection with Shift held. Ruled (a).
- **MQ4 — arcs.** (a) A fillet of radius r is n = max(2, ⌈θ/15°⌉) segments over its turn angle θ (6 for a right
  angle); a cylinder 16 sides; both on the 0.1 nm quantum; where a segment is shorter than the generator's 6 s
  trace, fewer segments and the reason shown. (b) A user-set count. **Recommended (a)** with the count editable in the
  inspector. Ruled (a).
- **MQ5 — Refine box.** Checked while proposing: both generators already take `SizeBox` lists
  (`MeshOptions2D::boxes`, `MeshOptions3D::boxes`), and per-region `region_max_size`, but the deck exposes neither.
  (a) A new opt-in `[mesh] refine` key mapped to `SizeBox`, so Refine box ships in slice 4 with no generator change.
  (b) Leave it out. **Recommended (a)**. Ruled (a).
- **MQ6 — the 2D editor.** (a) Gets MT1–MT3 and the Edges group's corner tools and the overlap rule; the 3D-only
  tools are hidden there. (b) 3D only. **Recommended (a)**: one way of working in both. Ruled (a).

**Risks:** the 3D generator's rules (26l SQ6) refuse some layouts. At the FinFET decks' s = 0.08 nm a 6-segment
fillet needs r ≥ 23 s ≈ 1.9 nm and a 16-sided cylinder r ≥ 15 s ≈ 1.2 nm (arithmetic, not measured), so sizes are
rarely the problem; MT12 keeps drawn layouts inside the rules. What stays refused: a nanowire through a gate (a holed
interface piece), fins or fillets touching the domain's sides, a plane crossing a slanted wall, sectors under 60°.
The first three are Unit 64; the last stays until measured (Unit 16 U12, the owner's ruling). Every non-right bend
adds graded pairs along its whole length, so fillets may cost many nodes (MT12's estimate shows it); the floating bar and live deck line touch every drag path, so MC4 and MC9
are the gates there; `[parameters]` is a deck change, so MC1 and MC7 guard every existing deck.

**Slice 0 as built (MT1–MT3, 2026-11-21; owner "start"; Release, MSVC 14.51, 6 cores).** `app::modeller`
(`include/NiTCAD/app/modeller.hpp`, no Win32): the groups and their tools, `navigation` (what a press does to the
view, camera free or fixed), `place_bar`, `deck_line`, `Entry` (typed values). `ModellerBar` (`src/app/modeller_bar.*`,
shared by both pages): the group switcher and the corner's buttons (Camera free/fixed (L), Reset view, Fit, Measure,
Mesh) in the strip where the fixed tool row was, the floating bar (the armed group's tools and its line) and the value
chip as children of the page over the canvas (the canvas clips siblings). Both pages: the camera lock on every drag,
the Pan tool, Space and the wheel; picks open their group; the status bar shows the deck line of the drag (the edit's
own writer on the in-progress model); the chip its value; typed digits (nm, um, cm; µm by default) set it, in the
drag's direction, Enter applying, Esc cancelling. Typed values are taken by the 3D axis, face, move, plane, base and
height drags and the 2D move and rectangle; on a two-axis drag they set the axis of the larger offset (recorded).
- MC1: the 26k and 26l window tests pass unchanged (789 assertions in 22 cases with slice 0's; no test named the old
  row); app unit tests 110/110. Decks untouched.
- MC2, MC3, MC4: `tests/app/modeller_window_test.cpp` (6 cases) and `tests/app/modeller_test.cpp` (5 cases, 6000
  random bar placements over 96/144/192 dpi). In the window table (3 dpi × 6 views) the bar never needed the corner.
- MC9: deck line + bar placement with 100 solids, median 0.197 ms (≤ 1 ms); SE5 orbit 5.10 ms, drag 5.21 ms; GE5 2.93 ms.
- MC11: 14/15 mutants caught; the survivor (the bar's overlap test removed) is equivalent, since a place beside the
  selection never overlaps it.
- Not done in slice 0 (recorded): Section stays in the inspector rather than the view's corner; Shift+left-drag pans in
  3D (MT1), the 2D editor keeps Shift for the rubber band (MQ3); the camera state lives in the page (kept across
  reloads and tab changes, not written to the view settings file). Whether the floating bar draws cleanly over the
  D3D12 flip-model canvas needs the owner's eyes (no screen automation).

**Slice 1 as built (MT4, MT5, 2026-11-21; owner "start slice 1"; Release, MSVC 14.51, 6 cores).** MT4: `[region]
overlap = replace|fill` (layered meshes only; elsewhere a build error), `LayeredShape::fill` / `LayeredSolid::fill`.
Both kernels paint in one order, replacing shapes top-first and then filling shapes bottom-first, the first to hold a
cell taking it (the 3D slanted labeller: the top replacing solid of a slab, else its lowest filling one). The views
draw and hit-test in the same order (`paint_order` in `app::geo` and `app::solid`). The 2D and 3D models read the
key, and write it by minimal edits: after the material of a new region; replaced, removed or inserted after the
material of a kept one. Set from the Solid group's Fill / Replace (fill unless every picked region already fills) or
the inspector's switch. MT5: the Transform group is Move, then the modes Copy, Mirror, Array, Align (`Mode`; the
armed one is `armed_key`). Copy (or Alt+drag in any mode) moves copies, the originals kept; a copy not moved, or
cancelled by Esc, is none. Mirror is a click: 3D about the centre across the clicked face's axis (Ctrl: across that
face's own plane; Shift: a quarter turn about the view's axis); 2D left-right (Ctrl: up-down; a side handle: across
that side; Shift: a quarter turn). Array is a drag from the solid: the pitch, the wheel setting the count (2-100,
whatever the camera), one region per copy keeping the material and `overlap`. Align is a click: the picked solids'
(or shapes') matching face goes to the clicked face's plane (a domain face too; 2D: the side of the clicked bounds
nearest the click). Model calls: `set_overlap`, `mirror_*_at`, `rotate_*`, `array_*`, `align_*` in both models. Each
is one deck edit and one undo step.
- MC5: `[overlap]` in the mesh tests: 1000 seeded layouts with fill shapes in 2D (seed 26131), in 3D rectilinear
  (26132) and slanted (26133), each against a sequential oracle, exact; with a fill shape over background only, the
  tiling is bit-identical to `replace`. Deck build: a fill box takes only what no earlier shape holds; the errors
  (second key, bad value, non-layered mesh). The mesh `[layered]` tests: 30/30.
- MC6: app tests `[overlap]` and `[transform]` (4 cases: read/write round trip, `only_geometry`, exact volumes, builds
  of arrayed, mirrored, turned and aligned layouts; four quarter turns give the model back; an L prism turned keeps
  its volume). Window tests (2 cases, by messages): Copy, Mirror, Align, Array with the wheel, Fill / Replace and the
  inspector switch in 3D and 2D, every step undone to the byte and redone.
- MC4: the deck line shown during the copy and array drags equals the line written (both window tests).
- MC1: the 26k, 26l and slice 0 window tests pass unchanged (910 assertions in 24 cases with slice 1's); app unit
  tests 114/114; deck tests pass. One test changed: `modeller_test`'s tool table, which names the group sizes (Solid
  4, Transform 5) and the Duplicate button (removed; Copy replaces it).
- MC9: deck line + bar with 100 solids median 0.159 ms (≤ 1 ms); SE5 orbit 4.81 ms, drag 4.85 ms.
- MC11: 16/16 mutants caught (kernel orders 2D and 3D, deck build, readers and writers, the resolver, Array's count
  and fill, Align's dedupe, Copy's cancel and unmoved press, the wheel during Array, Copy mode). Runs 4-16 sharded
  6 ways after a first sharded run was void (its runner never started the tests; recorded).
- Not done or limits (recorded): `overlap` is per region, not per solid; Alt+drag is not tested by messages (it
  reads `GetKeyState`); Align translates, it does not resize; the quarter turn rounds its centre down when it falls
  off the 0.1 nm quantum; the 3D view draws fill solids depth-tested, not in paint order; a mirror or turn leaving
  the domain is refused by the writer (the deck kept); the 2D model builds of the transform cases use a plain
  silicon deck, since small boxes near the layered deck's well corner trip the generator's P4 corner rule (U12).
  In the 3D view only top faces of solids under the oxide are pickable, so a buried solid is reached by the
  inspector or the selection.

**Slice 2 as built (MT6, MT7, MT12, 2026-11-22; owner "start slice 2"; Release, MSVC 14.51, 6 cores).** Shared corner
geometry in `app::geo` (`arc_segments`, `polygon_sides`, `corner_turn`, `fillet_tangent`, `chamfer_corner`,
`fillet_corner`, `regular_polygon`; 2D `chamfer_shape`, `fillet_shape`), on the quantum, a cut point landing on its
edge's far vertex merged with it. 3D model (`app::solid`): `edge_corner` (an edge along an axis as a prism's corner, a
box written as its rectangle), `chamfer_edge`, `fillet_edge`, `fillet_vertex`, `taper_edge`, `add_cylinder`,
`face_normal`, `face_solid`, `inset_face`, `extrude_face`. MT12: `mesh::check_rules` (generate_mesh's rules only, no
node placed; its body now `generate(…, rules)`, generate_mesh unchanged) and `mesh::largest_resolution` (down by
eighths from the complex's extent to a pass, then bisection to 0.1 %); `solid::automatic_resolution`, `keep_clear`,
`extent` (volume, interface area, feature-line length), `node_estimate`, `mesh_length`, `set_mesh_length`.
Pages: the Faces group (Push/Pull, Inset, Extrude, Contact), the Edges group (Chamfer, Fillet, Round vertex, Taper),
Cylinder (O) in Solid; `modeller::group_of`, `tools_for` (the 2D bar leaves out the 3D-only tools: MQ6 gives 2D
Chamfer and Round vertex on a shape's corners). Inset is two steps (drag the margin, release, move for the depth,
click); Cylinder three (centre, radius, height). The inspector's "Use the largest interface_resolution" writes the
automatic s (rounded down to the quantum) as one edit and reports the node estimate.
- MC6: app tests `[edges]`, `[faces]` (arcs and corners; Chamfer, Fillet, Round vertex, Taper, Cylinder and Extrude each
  written by minimal edits, read back, undone to the byte, the partition's volume exact (to 1e-12) and the deck
  building; Inset: see limits). Window tests `[tools]` (2 cases, 127 assertions): every tool in 3D and Chamfer and
  Round vertex in 2D by messages, each deck equal to the model function's own edit, every step undone and redone.
- MC4: the deck line shown during every slice 2 drag equals the line written (window tests).
- MT12 (2): fillet chords and cylinder sides drop while one is under 6 s (the chip says how many, or that they are under
  6 s); Copy and Array offsets keep clear (onto a plane or 4 s off; `keep_clear`); Inset's margin and depth and
  Extrude's depth are 0 or at least 4 s.
- The automatic s is `mesh::largest_meshable_resolution`: largest_resolution's s (the rules only), then built; a
  refusal a finer s can cure (degenerate_mesh: the interface certificate) steps s down by a fifth a build, at most 4
  builds (about half the rules' s; each finer build costs more), then a bisection on builds to 0.5 %; any other
  refusal (invalid_input: a pair set's geometry) ends it at once. Cost on the 7-layout table: 0.9-4.9 s a layout (one
  build each: the rules' s meshes there).
- MC13 table (7 drawn layouts): the automatic s meshes and s × 1.01 the generator refuses; the rules' own s passes
  check_rules and its 1.01 does not, each layout.
  Node estimate N = 2.65 V/h³ + 1.12 A/hp² + 13.3 L/hp (hp = min(h, 4 s)), fitted by least squares on that table's 7
  meshes (7k-15k nodes): within ±10 % there. In-sample only: no held-out table (recorded).
  500 random layouts (seed 26231, hidden test `[mc13]`, 6 processes by NITCAD_SHARD): accepted 301, all meshing (0
  false accepts); refused 81 (118 more stop in the partition: slanted prisms along two axes). **MC13 missed** on two
  counts, recorded, not loosened: 4 false refusals (the 4-build cap gives up where a finer s meshes) and 4 layouts
  that also mesh at 1.01 s (meshing is not monotone in s once the certificate counts). A first run with the rules-only
  s had 19 false accepts (P4b certificate) and a crash (layout 284): both fixed, below.
- Bugs found and fixed in slice 2's checks: (1) the rules-only s was refused by the certificate after meshing (19 of
  500): the build-confirmed search above. (2) Layout 284 (a cylinder whose axis lies on a box face) at s = 0.5 nm: a
  pair set's 2D mesh held two copies of one vertex an ulp apart, and splitting the segment between them inserted a
  point twice, ending in Triangulation2D's contract (`in_circle`, existing generate2d code): generate2d now refuses
  a point already a vertex (`coincides`, the first insertion, segment splits, circumcentres) as degenerate_mesh;
  regression test `[l284]`. (3) Inset's depth was measured from the release pixel's ray, not the face: after the
  margin drag, 0.1um straight up read 0.15um (view-dependent); now from the face (found by mutant 16). (4) At 200 %
  the Solid group's five buttons overflowed the canvas: the bar narrows its buttons to the viewport (48 px floor at
  96 dpi).
  (5) mesh_length read the first [mesh], set_mesh_length the last: both the last, as read_model. (6) Replacing a
  continued entry at the end of a deck without a final newline added one: kept.
- MC1: app tests 45 570 assertions (full app binary) and `[solid]` 10 744/22 after the fixes; window tests
  `[modeller]`, `[solid]`, `[geometry]` 1231/0; mesh `[complex_generate]`, `[complex_interfaces]`, `[layered]` all
  pass; deck tests pass (9 expected failures). clangd (`--check`) and clang-tidy (bugprone, analyzer) on the changed
  files: no error in slice 2's code (one raw CR in a test char literal fixed).
- MC9: slice 2's per-frame work with 100 solids: a Copy frame (keep_clear and its deck line) median 0.152 ms, a
  Fillet frame 0.148 ms (≤ 1 ms); slice 1's 0.149 ms; SE5 orbit 4.97 ms; GE5 2.67 ms.
- MC12: met (no Win32 in the models, no app or deck in mesh, no new dependency).
- MC11: 20/20 caught (~/mut-gv/s26m2.sh, s26m2b.sh through run6.sh, sharded on 6 cores). First pass 17 of 21: three
  did not compile (an unused variable under /WX; rewritten) and "a handle press keeps no start point" survived. It
  was no fault of the tests: a handle's pick already carries its place on the face, so the page code it mutated
  (added for a test that could not find a face under the handles) did nothing; removed, so that mutant is gone. The
  window tests now check the dragged value before anything is typed (Extrude, Inset's depth), which found bug (3).
- MC1 (after the last fixes): window `[modeller]`, `[solid]`, `[geometry]` 1233/27.
- Limits (recorded): Inset with a margin on every side of a face inside the device leaves that face a ring, a holed
  interface piece this branch's generator refuses (Unit 64, merged on main, is not on this branch); the refusal shows,
  and the automatic s says it. Inset of the face's own material joins its region. Pushing a solid's faces needs them
  pickable: solids inside the domain are hidden by its faces (26l), reached with the section plane. Align is not
  kept clear (it places a face exactly). `check_rules` leaves out what only a mesh shows (the pair sets, the
  certificate): the automatic s builds to confirm it, up to a few seconds a layout (more where it steps down). The
  automatic s is shown on demand (the button), not live in the chip; the node estimate shows in the button's status,
  not during drags. Arrays keep clear on the first copy's offset only. MC13's 500-layout test is hidden (it meshes
  every layout; about an hour on 6 processes).


**Slice 2 follow-up (2026-11-22, owner "fix").** The automatic s climbs after its search: while 1 % above still builds,
s rises by 1 % (bounded: above largest_resolution's s the rules refuse, without placing a node). MC13 rerun (500
layouts, 6 processes): accepted 301, all meshing (0 false accepts), **0 meshing at 1.01 s (was 4)**; false refusals
still 4 (layouts 124, 381, 463, 473: each of the 4 builds fails the certificate, P4b "the dual face of an edge in a
corner zone", where the test's 0.005um meshes) — **MC13 still missed on that count, recorded**: the cure is in the
generator's corner zones, not the search. The test's message comparison (a refusal by the rules at the automatic s
against the rules at 0.005um) is informational, not MC13: which rule fires first depends on s. The climb's mutant is
caught only by the hidden MC13 run (an hour), not run as a mutant. Two range-for copies in the window test fixed
(clang -Wrange-loop-construct).

**Slice 3 as built (MT8, 2026-11-22; owner "start slice 3"; Release, MSVC 14.51, 6 cores).** Deck: `[parameters]`
(once, no name) of `name = value`, a value a length (with nm, um or cm), a number, or a sum of lengths and other
parameters; references `$name` and sums `$a + $b`, `$a - 5nm`, `5nm + $a` in [mesh], [region], [contact] and [doping]
(MQ1 (a); [contact] added to MT8's list: a patch box is geometry too). An operator has spaces on both sides or on
neither: `$a -5nm` is two values. `deck::substituted` replaces them before the build reads values (deck::build calls it,
the plan keeping no view into the deck; [adapt] substitutes per pass): a lone reference by its parameter's text as
written (the same bits as the value written out), a sum by its value in cm written shortest round-trip. Errors on the
referring line (invalid_input): undefined, a cycle (`$a -> $b -> $a`), a number in a sum, a value neither a length nor
a number, a bad or repeated name, a second [parameters], a reference outside those sections, a `$` not starting one.
A number parameter in a length place reads as cm, exactly as its value written out would (DK5: a unitless length is
cm): the "non-length where a length is needed" of MC7 is a number in a sum and a value neither (recorded).
Modeller (`app::solid`): `Model::parameters` (name, expression, resolved), `Solid::refs` and `Model::domain_refs`
(each coordinate's reference as written); read_model resolves them; write_model writes a coordinate as its reference
while the reference gives it (else its value: the tie broken), [parameters] by minimal edits (a changed expression in
place, a new one after the last, a new section before the first). `set_parameter` (every tied face follows),
`tie_edits` (after an edit: a moved coordinate tied by a lone `$p` sets p; one tied by a sum moves its last
parameter by the rest; the edited solid's other tied coordinates are held where the edit put them, so dragging a box's
near face keeps its far face and Lg takes the rest; a tie to a parameter that is itself an expression, or that
disagrees, is broken), `name_extent` (N: the extent named, the far face `<near> + $name`, a low face `<near> - $name`).
A parameter a drag sets is written with a unit (a bare 0 would be a number). Page: the inspector's Parameters (an
editable `name = value` list and Apply parameters, one undo step; refused with the reason in the status when one
does not resolve); commit_model and the live line follow ties (MC4); N during a face drag, nothing typed, names the
dragged extent: a new parameter L<k>, then each existing length parameter in turn (the chip says which); once a number
is typed, n is its unit's.
- MC7: deck tests `[parameters]` 88 assertions: the diode (mesh, doping, number parameters) and a layered generate2d
  deck (domain, solids, sums, tied faces) build to the same run identity and mesh digest as the values written out;
  every error on its line. App `[parameters]` (model): ties read and written back byte for byte; Lg changed moves both
  tied faces with one line changed; errors leave the model as it was; drags of the far face, the near face, the whole
  solid, a low face named with `-`; an expression parameter's tie broken; N naming. Window `[parameters]` 41: Lg in the
  inspector moves every tied face (one undo step, a cycle refused), a tied face dragged by typing (the line shown the
  line written, MC4), N naming and cycling.
- MC9: a tied drag frame with 100 tied solids (the edit, tie_edits, write_model, the deck line) median 0.68 ms
  (≤ 1 ms; measured while MC13 ran on 3 cores); window suites' MC9 0.19 ms, slice 2's frames 0.22/0.18 ms.
- MC1: app tests 45 755 assertions in 122 cases (before the last tests were added; `[parameters]` 387 after); window
  `[modeller]`, `[solid]`, `[geometry]` 1271 in 28 (then `[parameters]` 41); deck 128 cases, 123 passed and 5 failed as
  expected (the 9 `[!shouldfail]` assertions, as before).
- MC11: 17/17 caught (~/mut-gv/s26m3.sh, s26m3b.sh through run6.sh, sharded on 6 cores). First pass 15: "a sum of
  lengths not evaluated" survived (only the app tests held one; a deck test added) and one did not compile (an unused
  variable under /WX; rewritten).
- Deck `[parameters]` 93 assertions after it.
- Limits (recorded): a contact's patch box reads references but is written as values (its ties are not kept); 2D
  (the geometry editor) does not read references yet; renaming a parameter in the
  inspector makes a new one (the old stays, its references with it); the live preview moves only the dragged solid,
  other tied faces following when the drag ends (the deck line shows the parameter's line throughout); a polygon
  vertex inserted or removed drops that prism's ties (its slots renumbered).

**Slice 4 as built (MT9, MT10, 2026-11-22; owner "once slice 3 is done start slice 4"; Release, MSVC 14.51, 6
cores).** Deck: the `[doping]` profile reader moved out of the builder into `deck::parse_profile` (the build reads
every profile through it, unchanged: deck baseline 123 passed, 5 failed as expected, the 9 `[!shouldfail]` as
before), with its writer `deck::profile_text` (lengths in cm, shortest round-trip: read back to the same doubles;
windows only where bounded). Modeller (`app::solid`): `profiles(deck)` (references substituted), `add_profile` (one
`profile =` line last in [doping], a new [doping] before the first contact or analysis section when there is none;
every other line byte for byte), `region_profile` (constant over the solid's bounds; the deck windows two axes, so the
unbounded depth axis is the one the solid spans most of the domain along, z on a tie: exact where it spans all of
it), `implant_profile` (gaussian along the face's axis, its peak `depth` inside, σ half of it), `junction_profile`
(erfc from a face looking down an axis, its length half the depth and its centre so the profile equals the other
dopant's net `depth` inside the face, at the window's centre: c = x_j - w erfc⁻¹(2 N_bg / N_peak), erfc⁻¹ by Newton to
the last bits), `profile_extent` (as drawn), `size_boxes` / `add_size_box`. Page: the Doping group's Region profile
(click a solid), Implant (drag into a face: the peak depth, typed or dragged), Junction (drag into a face looking down
an axis: the junction depth); D swaps donor and acceptor and the wheel sets the level by decades (1e12 to 1e22),
during a drag too; the chip shows the depth and the dopant and level. The Mesh group's Max size (drag up or down
anywhere: doubled every 100 px, or typed; the chip estimates the nodes, MT12) and Refine box (drag a rectangle on the
work plane: a box through the domain along the plane's normal, its size max_size / 4, the wheel halving or doubling
it). Each writes one line by minimal edit, shown live as the deck line (MC4), one undo step. Profiles are drawn as
see-through boxes (donors blue, acceptors red, deeper with the level, log 1e14 to 1e21) while Doping is open, refine
boxes grey while Mesh is open, the dragged one too; both read once per deck change, not per frame. Snapping (found
building MC10): a snap onto a face or an edge now also puts its free coordinates on the view's grid where that stays
on the face or edge (it kept the pixel's place: a gate at 0.0401um beside a fin at 0.04um, slivers of a quantum the
generator refuses).
- **MQ5 premise corrected:** the deck already exposes the generators' size boxes, `[mesh] box = u0 v0 u1 v1 h`
  (generate2d) and `box = x0 y0 z0 x1 y1 z1 h` (generate3d, layered too). Refine box writes that existing key; no
  `refine` key was added (recommended option: no second spelling of one thing).
- MC8: app `[doping]` 3049 assertions: each tool's line exactly (a minimal edit: one line more, the rest byte for
  byte), read back to the same profile, the junction meeting 1e16 to 1e-9 at its depth, its three refusals; max_size
  and a refine box; then the deck built and the net doping at every one of its 1506 semiconductor nodes equal to the
  profiles' formula (`device::evaluate` summed) to 1e-12 relative. Window `[doping]` 61: Region profile (wheel, D, a
  click), Implant (typed), a Junction from a face looking up refused with the deck as it was, Max size (the chip's node
  estimate), Refine box (the wheel's size); each the line shown, one undo step and redo; the deck builds.
- MC10: a FinFET by window messages only (`[finfet]`): Box (a silicon fin), Array (two, pitch typed), Fillet (each
  fin's top edge, 0.02um), a sio2 gate box with Fill and its x_max face dragged with N naming `Lg` (`0.04um + $Lg`),
  a Gaussian implant into each fin's top (donors 1e20, 10nm), Region profile (acceptors 1e17 in each fin), three
  Contacts (source and drain the x faces, the gate on top made an electrode by the inspector's Kind), the inspector's
  largest interface_resolution (1.2nm); built, **13 359 nodes** (≤ 20k), equilibrium converged in 1.7 s.
- MC9: 100 solids and 100 profiles: an Implant drag's frame median 0.12 ms, a Junction's 0.11 ms, the 100 volumes
  drawn +0.23 ms (≤ 1 ms).
- MC1: app tests 49 059 assertions in 124 cases (the slice 2 group test updated: Doping has its tools, Mesh lacks
  Interface layers); window `[modeller]`, `[solid]`, `[geometry]` 1386 in 30, all passing after the snap change.
- MC11: 12/12 caught (~/mut-gv/s26m4.sh through run6.sh, sharded on 6 cores): windows not read, lateral not
  written, the region depth axis always z, the implant peak outside, the junction centre off, a junction from a face
  looking up, a new [doping] section every time (its first pattern did not apply; rerun), a refine box unchecked, the
  drag's line not written, D doing nothing, the wheel keeping the level, face snaps off the grid.
- **Misses and deviations (recorded):** Interface layers (MT10) **not built**: the `layer` key exists only for
  extrude3d cross-sections (`option_swept`); a layered deck has none and MT10 allows no generator change — the Mesh
  group says so. Implant's window is the whole face and its lateral spread 0 (MT9: a window dragged on the face, then
  the depth and the spread). Region profile's level is the wheel (MT9: a drag up or down), the same for all three
  tools. Refine box spans the domain along the work plane's normal and takes its size from the wheel (MT10: a box
  dragged, then inside it for the size). MC10's implant covers each fin's whole top, under the gate too (no mask: the
  window is the face), so the device is a resistor-like FinFET; MC10 asks only that it builds and runs equilibrium.
- Limits: a junction is made only from a face looking down an axis (the deck's erfc falls along its axis); the
  junction's background is taken at the window's centre only; profiles are drawn only while the Doping group is open;
  the profile tools are 3D only (the 2D editor has no Doping group tools yet, MQ6 allows it).

**Slice 4 follow-up (2026-11-22, owner "do the follow ups").** Three of the recorded deviations built as MT9 and
MT10 give them: Implant is a window dragged on the face (snapped on the face's plane as on the work plane, kept on the
face; a click without a drag takes the whole face), released, then the peak depth (moved or typed) and a click, then
the lateral spread (the distance out from the window on the face, or typed) and a click writing the line. Region
profile is a press on a solid dragged up or down for the level, on the 1, 2, 5 steps (3 a decade, one every 30 px;
a click keeps the level; the wheel still steps decades and D swaps the dopant). Refine box is a rectangle dragged on
the work plane, released, then its height along the plane's normal (moved or typed) and a click, then its size (up
doubles, down halves, or typed; the wheel too) and a click writing the `[mesh] box` line. Enter writes any phase's
line as it stands. A press in a Doping tool picks a solid's face under the handles (a picked solid's gizmo moved it
before: found by the new window test) and never falls back to moving (a press on empty background still navigates, MT1). MC10 now masks the gate: two windows on each
fin's top (x 0..0.04um and 0.08..0.12um), the peak 10nm in, lateral 2nm: source and drain only, six profiles; built,
13 359 nodes, equilibrium 2.0 s.
- Window `[doping]` 135 assertions in 2 cases (an edge press in Implant moving nothing; each phase's chip, the deck unchanged until the last click, the line
  shown the line written, one undo step); app 49 059 in 124; window `[modeller]`, `[solid]`, `[geometry]` 1604 in 31.
- MC11: 7/7 caught (~/mut-gv/s26m4b.sh, three rerun by hand with fixed patterns): the level drag doing nothing, the
  window always the face, the lateral spread not applied, the depth click ending the implant, the refine height
  ignored, the window off the face's grid, a press on an edge falling back to moving (survived first: the edge-press
  test added).
- Still a miss: Interface layers (no `layer` key for a layered deck without a generator change).

**PR #116's first CI run (2026-11-22): two failures, fixed.** MC3: at 200 % the bar was wider than the view on
the runner's 1024 x 768 screen (a 525 px view, the bar 602 px: the hint's 300 px at 96 dpi a floor). The hint now
narrows to the view's room and wraps onto as many lines as it takes, and the buttons wrap into rows when even their
48 px floor leaves no room. MC3 now runs at three window sizes (as opened, 1024 x 768, 800 x 600: 48 placements; the
old bar fails at 1024 x 768 exactly as on CI); 800 x 600 at 200 % is left out (a 301 x 195 px view holds no bar in
any layout: a limit). The test's check of the bar's place is now to 1e-9 px (the page places in its panel's
coordinates; the new sizes reached "above" placements, differing by rounding). MC9 slice 3's tied drag frame took
1.055 ms on the runner: recorded there, gated on this workstation (0.769 ms), as SE5 (owner, 2026-11-22). Window
`[modeller]`, `[solid]`, `[geometry]` 1722 in 31; app 49 058 in 124. Seen once: MC2's L key ignored in one full run
while the owner typed (the page reads Ctrl from the physical keyboard, GetKeyState); it passed in the next full run and
3 of 3 alone: a test hazard, not changed.

### Unit 26n — devices from dimensions: the dimensions panel, a parametric device library, dimension sweeps (proposed 2026-11-23; on `app/device-templates`)

**Owner request (2026-11-23):** "Edit in Geometry is difficult … make the Geometry Section easy for creating
semiconductor devices and research purpose"; then "investigate NiTCAD's existing Geometry system and propose the best
solution to make semiconductor device creation simple, intuitive and research-oriented"; of the three approaches put
(A devices from named dimensions with templates, B a process-style layer stack, C polishing the freehand modeller
only) the owner chose **A** ("go with A").

**Where things stand (found while proposing).**
- The 3D modeller (26l, 26m) and the 2D editor (26k) are general CAD: drag boxes, prisms and cylinders on a work plane,
  then Faces/Edges/Transform tools. Nothing knows a gate, a fin, an oxide or a source.
- A device's numbers are written out by hand, many times over. In `examples/mosfet_3d.deck` the gate edge 0.5 µm is
  in the gate contact's box, the source implant's window and (implicitly) the mesh lines; the oxide's 0.95–1 µm is in
  the region, three contacts and the mesh. Changing Lg is about five edits that must agree.
- `[parameters]` with `$name` references (26m MT8, `src/deck/parameters.cpp`) already substitutes in `[mesh]`,
  `[region]`, `[contact]` and `[doping]`, but no shipped deck uses it and the only editor is a bare text box in the 3D
  inspector ("Parameters (name = value; $name in solids)").
- The modeller edits only `generate3d` layered decks (`MeshKind::layered`); every shipped deck is `tensor` (7) or
  `extrude3d` (4), so the shipped devices cannot be edited in it at all, and 1D/2D tensor decks get no editor.
- `[sweep]` steps contact biases only. Comparing Lg = 20, 30, 45 nm means three hand-made decks, and the comparison
  view (26b) holds two runs (A and B). Unit 24's figures (V_th, SS, g_m max, on/off) are extracted per run.
- The tester's screenshot of 2026-11-22 (cylinders drawn outside a 722-node domain, "region has no nodes") is the
  freehand path's typical failure: nothing keeps a drawn solid inside the domain.

**Approaches considered:** (A, chosen) a dimensions panel over `[parameters]` plus a library of parametric decks,
with sweeps over a dimension; (B) a process-style stack (substrate, deposit, etch, implant) — the most general, but a
new subsystem imitating a process simulator NiTCAD does not have, and 26m's owner ruling left process emulation out;
(C) polish the freehand modeller only — needed anyway (slice 4 here), but it does not remove the repeated numbers or
give sweeps.

**Standing conditions:** as Units 26l/26m. The deck is the only source of truth, written by minimal edits through the
deck editor's undo. Every deck that builds today builds to the same bits. Win32 + D3D12 + Direct2D only, no new
dependency, no Python. The models have no Win32 type. Criteria fixed before code. Targeted Release tests per slice;
full suites and Debug only on the owner's word; builds `-j 6`, runs on 6 threads. No screen automation. Nothing
pushed or merged without authorization.

**Left out:** process emulation (owner, 2026-11-21); a constraint solver (dimensions are independent parameters and
sums of them, MQ1 of 26m); GDS/STEP import; optimisation loops (sweeps only); remote or parallel batch runs (one run at
a time on 6 threads); new physics.

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| DT1 | 0 | **The dimensions panel.** For any deck with `[parameters]` (1D, 2D, 3D; tensor, extrude3d, generate2d/3d), the Geometry tab's inspector shows a table in place of the text box: one row per parameter — its label, value with unit, and the kind (length / number). A label is the parameter line's trailing comment (`Lg = 20nm  # Gate length`), else the name. Typing a value and Enter (or the mouse wheel over the cell: ×/÷ 10^(1/10) for numbers, ±1 lattice step × a 1-2-5 ladder for lengths) writes the one `[parameters]` line by a minimal edit; the live check (`live_check`) runs and the preview rebuilds. A refused value shows the build's message on that row and is not written. Rows may be grouped by a comment line `# group: Gate` inside `[parameters]`. | Decks without `[parameters]` show "No named dimensions — Name this value…" and the 26m tools are unchanged. |
| DT2 | 0 | **The Geometry tab for every deck kind.** A deck the modeller cannot edit (tensor, 1D, 2D tensor) shows the device preview (the existing build preview, "Device (preview, not solved)") with the dimensions panel beside it, in place of today's blank or Device-view-only state. 1D shows the layer bar with dimension labels; 2D/3D the existing views, each dimension drawn as a labelled dimension line between the faces it sets (a parameter's faces are found by substituting it ± one lattice step and diffing the built geometry). | Clicking a dimension line selects its row; dragging it changes the value (one edit on release). |
| DT3 | 1 | **The device library.** Parametric decks in `examples/` (each a full runnable deck: geometry, doping, contacts, mesh, models, a sweep, outputs) for: PN diode (1D, 2D), Schottky diode (1D), MOS capacitor (1D, C–V), planar MOSFET (2D, 3D), SOI MOSFET (2D), FinFET (3D, bulk and SOI), JFET (2D). Every geometric number, doping level and junction position is a parameter; the mesh lines are written from the same parameters (`linspace($x0, $x0 + $Lg, n)`), so an interface always has mesh lines on it. Existing example decks are rewritten to their parametric form only where they build to the same bits (MC1); otherwise the parametric deck is a new file beside them. | 1D/2D tensor meshes keep sweeps cheap and exact; 3D FinFETs use `extrude3d` as today. |
| DT4 | 1 | **New device.** File ▸ New device… and the Start page's first card open a small dialog: the family list (with its dimension), a thumbnail, and the main dimensions (the panel's first group) to fill in before Create; Create copies the deck into the user's chosen folder (`copy_example`) with those values and opens it on the Geometry tab. | Same examples list the Start page reads (26 VQ11). |
| DT5 | 2 | **Unlock shapes (3D).** A button on the panel for an `extrude3d` or `tensor` 3D deck: converts it to a `generate3d` layered deck (26l SQ5's `convert_to_layered` for extrude3d; for a tensor deck its region boxes become box solids, its mesh lines' spacing becomes `max_size`/`interface_resolution` as the 26m MT12 estimate needs), keeping `[parameters]` references in the solids it writes, then opens the 3D modeller. A refusal says why (e.g. a layout the 3D generator refuses, Unit 64). | One-way and undoable (deck undo). |
| DT6 | 3 | **Dimension sweeps.** Each panel row has Sweep…: a value list (`20nm, 30nm, 45nm` or `from : to : step`, lengths or numbers). Run executes the deck once per value with that parameter overridden (no copy of the deck is edited; the override is applied after reading, as substitution is), writing each run's outputs to `<deck stem>_sweep_<name>/<name>=<value>/`. The same from the CLI: `nitcad run <deck> --set Lg=30nm` (repeatable) and `nitcad sweep <deck> Lg=20nm,30nm,45nm`. | One run at a time; Stop ends after the current point. |
| DT7 | 3 | **Sweep results.** The Results tab gains a sweep view: the runs' curves overlaid (one colour per value, viridis order, legend `Lg = 20 nm`…), and Unit 24's figures plotted against the swept parameter (V_th, SS, g_m max, on/off vs Lg) with a CSV export `<…>_sweep_<name>/figures.csv`. Run comparison (26b) still takes any two runs. | No new figure extraction. |
| DT8 | 4 | **Modeller guardrails.** Drawing in the 3D modeller and the 2D editor keeps solids inside the domain (a drag clamps at the domain's faces; a solid moved out is refused with the reason, as a meshability refusal); the inspector's position and size fields for a picked solid take typed values and `$name` references; "Name this value" (26m's N during a drag) is also on the inspector field's context menu. | Fixes the 2026-11-22 screenshot's case. |
| DT9 | — | **Commits:** this proposal, a commit per slice, then ARCHITECTURE (6.14, the deck chapter for `--set` and labels, 11.3) and the records. The owner may stop after any slice. | |

**Criteria (fixed before code; a target missed is recorded as a limitation, never loosened):**

| # | Criterion |
|---|---|
| DC1 | **Nothing else moves:** every example and test deck not rewritten builds to the same bits; a rewritten example builds to the same bits as before (else it is not rewritten, DT3); the 26k/26l/26m suites pass unchanged except where a test names the parameters text box (recorded). |
| DC2 | **Panel edits are minimal:** for each library deck and each parameter, a typed value changes exactly that `[parameters]` line (byte diff of one line); undo restores the deck byte-identical. |
| DC3 | **One number moves everything:** for each library deck, for each length parameter p, the deck built with p + δ (δ one mesh-independent step) equals the same deck with every dependent value written out by hand at p + δ (geometry, doping windows, contacts, mesh lines) — compared on the built device's regions, contacts and net doping at 200 seeded points to 1e-12 relative. |
| DC4 | **The library runs:** every library deck at its default values builds, meshes under 20k nodes (3D) / 5k (2D) and runs its sweep to completion in Release on 6 threads within the existing slow-test budget; each MOSFET/FinFET's extracted V_th and SS are finite and SS ≥ 60 mV/dec × (1 − 1 %). |
| DC5 | **Sweeps are the runs they stand for:** `nitcad sweep d Lg=a,b` writes, for each value, outputs byte-identical to `nitcad run` on a copy of d with `Lg = a` written in; `--set` on an unknown name or a wrong kind is a usage error naming the parameter. |
| DC6 | **Dimension lines:** for each library deck, each length parameter's dimension line joins the faces the parameter moves (DT2's diff) and its label shows the parameter's current value; dragging it by n lattice steps writes the value + n steps. |
| DC7 | **Guardrails:** 500 seeded random drags of the Box, Prism and Cylinder tools and of the move gizmo never leave a solid outside the domain; a typed out-of-domain size is refused with the reason and the deck unchanged. |
| DC8 | **Unlock:** each 3D library deck either converts and builds as a layered deck whose regions at 200 seeded points equal the original's, or is refused with the 3D generator's reason (recorded per deck). |
| DC9 | **Speed:** a panel edit to the refreshed preview ≤ 300 ms median for the 2D library decks and ≤ 1 s for the 3D ones on WARP (Release); the modeller's 26l SE5 / 26m MC9 frame budgets unchanged. |
| DC10 | **Interaction:** by window messages only — New device ▸ FinFET (3D), set Lg = 30 nm and Wfin = 8 nm in the panel, Sweep… Lg = 20, 30 nm, run: two runs, the overlay shows two curves, figures.csv has two rows. |
| DC11 | **Mutants:** at least 8 per slice caught (e.g. a panel edit touching a second line; an override applied to the wrong parameter; a sweep folder reused between values; a clamp letting a solid one lattice step out). |
| DC12 | **Layering:** models in `app` without Win32 headers; overrides and labels in `deck`; the sweep driver in `deck`/`cli` without `app`; no new dependency. |

**Questions (each with a recommendation; ruled by the standing "choose recommended" rule, 2026-11-11):**

- **NQ1 — where a parameter's label lives.** (a) The parameter line's trailing comment, groups by `# group:` lines —
  no grammar change, readable in the deck. (b) New keys (`Lg.label = …`). **Recommended (a).** Ruled (a).
- **NQ2 — the library's mesh in 3D.** (a) Keep `extrude3d` (as the shipped FinFETs), Unlock converts to layered (DT5).
  (b) Write the library as `generate3d` layered so the modeller edits it at once — but at 1 nm oxides the 3D
  generator needs s ≈ 0.08 nm and refuses some fin layouts (Unit 64), which would make sweeps fail part-way.
  **Recommended (a).** Ruled (a).
- **NQ3 — arithmetic in parameters.** (a) Keep 26m MQ1's `+`/`−` of lengths; library decks are written so every
  position is a sum (e.g. `x_gate1 = $Lsd + $Lg`), no halves. (b) Add `×`/`÷` by numbers now. **Recommended (a)**;
  (b) if a library deck turns out to need it, recorded. Ruled (a).
- **NQ4 — where sweep outputs go.** (a) `<deck stem>_sweep_<name>/<name>=<value>/` beside the deck, each a normal
  run folder (opens in the results browser). (b) One combined file. **Recommended (a)**: every run stays an ordinary
  run, comparable with 26b. Ruled (a).
- **NQ5 — sweeping two dimensions.** (a) One parameter per sweep in this unit; a grid of two later. (b) Grids now.
  **Recommended (a).** Ruled (a).

**Risks:** DT3 is the largest slice (ten decks, each validated); a family whose physics needs a model not yet in
NiTCAD (e.g. a JFET's gate junction is fine; a Schottky barrier uses the existing `kind = schottky`) is dropped and
recorded rather than approximated. DT2's dimension lines rely on diffing built geometry, which is exact on tensor and
layered meshes but may find no face for a parameter used only in doping (then no line, the row only). DT5's tensor →
layered conversion may be refused for thin-oxide FinFETs by the 3D generator's rules; DC8 records it per deck. Sweeps
run whole decks serially, so a 3D sweep of five values costs five runs (the 18-min FinFET is not a library default:
DC4's node limits keep defaults quick).

**Slice 0 record (2026-11-23, branch `app/device-templates`, commits ee43bc8..HEAD, local only):** DT1 and DT2 as far
as slice 0 goes.

- **Built.** `app/dimensions` (rows from `[parameters]`, one-line edits keeping key and comment, the wheel's 1-2-5
  ladder, spans found by substituting at p and p + δ and diffing boxes, doping windows and linspace pieces, marks,
  `substituted_text`); text labels in field scenes (renderer D2D overlay); the dimensions panel in the Geometry
  inspector (Enter commits, Esc restores, wheel steps, refusals shown red under the row, disabled with the message
  while the deck does not parse). The 2D editor and the 3D modeller host the panel and draw the lines (2D: 12 px off
  the extent, before the rulers; 3D: on-top line batches with labels); a line is selected and dragged to change its
  value (one undo). The 3D inspector's parameters text box is replaced by the panel. Decks neither editor draws (1D,
  3D tensor) show the device preview in the canvas's place with the panel beside it; a 3D preview draws the lines.
- **Measured.** DC9 on this machine (Release): 2D panel edit → redrawn canvas median 2 ms (≤ 300 ms); 3D tensor
  panel edit → refreshed preview median 433 ms (≤ 1 s). DC11: 10 of 11 mutants killed (two after adding tests); the
  survivor (the panel's own unparsed-deck guard) is equivalent — `set_dimension` refuses the same decks. DC12: the
  model has no Win32 headers. DC1: `nitcad_deck_test "[deck]~[perf]"` 14468/14469 — the one failure is the hidden 26l benchmark "reduced FinFET as layered generate3d", which searches LF text in `decks/finfet_3d_reduced.deck` and this Windows checkout has it CRLF (`core.autocrlf = true`); the branch changes no deck, example or deck test. Targeted suites: `[dimensions]`, `[solid]~[perf]`,
  `[geometry]~[perf]`, `[modeller]~[perf]` pass; full Release/Debug suites not run (owner's word needed).
- **Rulings.** A bare number typed for a length takes the row's unit, so `1e20` on a length is `1e20nm` (the refusal
  test uses `20 nm nm`). A stepped number is written in shortest scientific form. A span is the shortest candidate
  with exactly one end moving; ties: boxes, then doping windows, then mesh pieces. A parametric 2D deck (shapes with
  `$` references) is shown view only from its substituted text: lines, drag and panel work, shape tools refuse with a
  status message. A selection's handles win over dimension lines in 3D; in 2D lines sit beyond the handles' reach.
  The line/label scene code is shared (`solid::add_dimension_marks`).
- **Limitations (recorded).** 1D decks show the panel but no lines (no layer bar yet). Polygons and prisms give no
  spans. A preview line selected in the panel re-highlights at the next refresh only; preview lines are not draggable.
  A dimension used only in doping gets a line only through its `x=a:b` window. A line on a continued (`\`) deck line
  is not substituted (view only refused, the editor's message stands).

**Slice 1 record (2026-11-23, branch `app/device-templates`, commits 083db06..HEAD, local only):** DT3 and DT4.

- **Built.** Ten parametric decks in `examples/devices/`: PN diode 1D and 2D, Schottky diode 1D, MOS capacitor 1D
  (C–V), planar MOSFET 2D and 3D, SOI MOSFET 2D, FinFET bulk and SOI (3D), JFET 2D. The MOSFET, JFET and SOI FinFET
  follow the owner's deck set (`decks/example-decks`, f839bab); the bulk FinFET follows the owner's
  `decks/finfet_bulk_3d.deck` (read only, not committed); the rest the shipped examples. Every geometric number,
  doping level and junction position is a `[parameters]` row, labelled by its comment, the first group `Main`.
  Positions are derived rows (`x_gate = $Lsd + $Lg`), and the mesh lines are written from the same rows. The four
  shipped examples are unchanged (the 26k tests edit `mosfet_2d.deck` in the 2D editor). File ▸ New device…
  (Ctrl+Shift+N, and the Start page's first button) opens a dialog over the window: the families, the chosen deck's
  thumbnail (its region boxes, a 3D deck cut at mid z), its main dimensions in a dimensions panel over the deck's
  text in memory (a value the build refuses shows on its row and is not written), and the folder (default
  Documents\NiTCAD\devices, Browse… picks one) and name. Create writes `<name>.deck`, never over a file, with its
  `[output]` files renamed after it, and opens it on the Geometry tab.
- **Measured.** DC2: every length row of every deck moved by one panel edit changes one line and still builds
  with the same node count. DC3: every geometric and doping number is a reference (a lint over
  `[mesh]`/`[region]`/`[contact]`/`[doping]`), and the 2D MOSFET at Lg + 1 nm equals the deck written out by hand
  (nodes, regions, doping, contacts, to 1e-12). DC4 (Release, solver threads 1): nodes jfet 4180, mosfet_2d 4331,
  pn_2d 2511, soi_2d 2257, bulk FinFET 6720, SOI FinFET 19343, mosfet_3d 8525; every deck runs at its defaults.
  mosfet_2d 11.7 s, V_th 0.446 V, SS 79.3 mV/dec; soi_mosfet_2d 3.8 s, 0.425 V, 92.8; mosfet_3d 57.2 s, 0.394 V, 79.6;
  bulk FinFET 25.8 s, 0.286 V, 84.2; SOI FinFET 1167.6 s, 0.229 V, 64.9; jfet_2d 15.2 s; the 1D decks under 0.2 s.
  **DC4 missed in part:** the SOI FinFET's run (19.5 min) is far over the slow-test budget. It stays a hidden test
  (`[.library_runs3d]`), not in ctest. DC11: 10 of 10 mutants killed (one after a new check). DC12: the model
  (`app/new_device`) has no Win32 headers.
- **Rulings.** The library is a folder of its own, beside the four examples, not a rewrite of them. The 3D FinFETs
  are tensor decks like the owner's, not `extrude3d`: the slice 0 preview draws a 3D tensor deck with its lines.
  The SOI FinFET's mesh along the fin is coarser than the owner's (17 lines, not 23) to stay under 20k nodes.
  Non-geometric keys (materials, work functions, Richardson constants, straddle widths) stay literals. DC4's
  figures are checked on the MOSFETs and FinFETs only; the JFET's ohmic gates give no transfer figures. The panel's
  Host gained a refusal hook (default: take the edit). The dialog closes from a posted message.
- **Limitations (recorded).** The dialog has no dialog keys yet (Tab, Enter and Esc are not routed through
  IsDialogMessage). Browse… is not driven by the tests. A thumbnail shows region boxes only (polygons and mesh
  files: the background). DC10's full interaction (with Sweep…) waits for the sweep slice.

**Slice 3 record (2026-11-23, branch `app/device-templates`, commits 1029021..HEAD, local only):** DT6 and DT7, taken
before slice 2 (Unlock) on the owner's word ("start the next slice (sweeps)").

- **Built.** `deck/dimension_sweep` (no `app`): `parse_override`, `with_overrides` (the parameter's line rewritten,
  key, spacing, comment and CR kept, the text parsed again with the deck's source and directory; a wrong kind is
  refused naming the parameter), `sweep_values` (a list kept as written, or `first : last : step` in first's unit),
  the sweep's folders and `with_outputs_in`. The CLI: `check`/`run <deck> --set name=value` (repeatable) and
  `nitcad sweep <deck> name=values [--set ...]`, every value checked before the first run. `RunController::start_sweep`
  runs the values one after another into their folders. In the app, each dimension row not derived from others has
  Sweep…; its dialog reads the values (a bare number taking the row's unit), says how many runs and where, and Run
  starts the sweep on the editor's text (an untitled deck saved first, as F5). When it is done the runs load
  (`app/sweep_results`), `figures.csv` is written beside them and the Curves tab shows "sweep: curves" (one curve
  per run, viridis in sweep order, legend `Lg = 20 nm`) and "sweep: <figure> vs <name>" for each of Unit 24's
  figures. `figure_curve` now holds the figures' choice of curve, shared by `extracted_figures` and the overlay.
- **Measured.** DC5: `nitcad sweep pn Lp=1um,1.5um` writes, per value, a device file and table byte-identical to
  `nitcad run` on a copy with the value written in, and a sweep file identical once its provenance is cleared (times
  differ by nature); `--set` likewise; an unknown name or a wrong kind exits 64 naming the parameter, and no folder is
  made. The GUI sweep's sweep files equal the CLI's to the bit (fields, currents, record). DC10 (hidden
  `[.sweep_dc10]`, Release, WARP): New device ▸ FinFET (bulk, 3D), Lg = 30 nm and Wfin = 8 nm in the panel, Sweep…
  Lg = 20, 30 nm, Run: two runs, two overlaid curves, `figures.csv` with two rows; 43.5 s. DC11: 12 of 12 mutants
  killed (two rewritten after their first form did not compile). DC12: the override and sweep code is in `deck`
  and `cli`; `app/sweep_results` has no Win32 header. Targeted suites pass: `[dimension_sweep]`, `[sweep_results]`,
  `[sweep_run]`, `[sweep]`, the run controller's, the dimensions, New device, window, geometry, modeller and solid
  window files, `run_test`, `parameters_test`, `syntax_test`, `workspace_test`, `compare_test`, `view_settings_test`.
  Full Release/Debug suites not run (owner's word needed).
- **Rulings.** Stop stops the current value's run at the solver's next safe point (its completed part written) and
  runs no further value, rather than waiting out the whole value (a 3D value can take many minutes); any run that
  stops early ends the sweep, as `nitcad sweep` exits 2. A sweep takes literal values only, and refuses a value
  given twice (two values would share a folder). A derived row has no Sweep… (it follows others). On the CLI a length
  needs its unit; in the dialog a bare number takes the row's. Any `--set` refusal is a usage error (64); a deck whose
  own parameters do not resolve fails as before (1). Lengths in `figures.csv` and on the figure plots are in nm.
  The overlay draws the curve the figures read (a gate sweep: the contact carrying most current; else the swept
  contact's own), on a log axis when the Curves tab's is. Opening a deck or a plain run clears the sweep views.
- **Limitations (recorded).** `nitcad sweep` writes no `figures.csv` (the figures are extracted in `app`; DC12 keeps
  the CLI's sweep driver out of it). The overlay takes sweep files only: a small-signal deck (the MOS capacitor's
  C–V) or a transient loads but shows no overlay curve. The sweep views are rebuilt from the folders only right after
  a sweep in the window; reopening a deck does not find an earlier sweep. The dialog has no dialog keys (as New
  device's).

**Slice 2 record (2026-11-23, branch `app/device-templates`, commits 7ba6205..HEAD, local only):** DT5, after slice 3.

- **Built.** `app/unlock` (no Win32): `kind_of` (a 3D tensor deck; an extrude3d deck) and `unlock`. A 3D tensor deck is
  rewritten by minimal edits: `[mesh]` becomes `type = generate3d`, `layered = on`, a `domain` from the mesh lines'
  ends as written (`0 $L 0 $y_ox 0 $W`), a `patch <contact> = ...` per contact on its face (the contact's box, its
  references kept; a contact without a box takes the whole face), `max_size` the widest mesh-line spacing and
  `interface_resolution` the largest the 3D generator's rules allow (`mesh::largest_resolution`, two significant
  figures, rounded down); the x, y, z and straddle lines go. Region boxes stay box solids byte for byte (one reaching
  past the domain cut to it, one wholly outside dropped, both in the note). Each contact's `patch = <face>` becomes
  `patch = <name>` and its box line goes. Refused, with the reason: an `[adapt]` section, a contact on a patch that is
  not a face or with two boxes, a layout the 3D generator refuses (its message), an s under 0.1 nm. An extrude3d deck
  goes through 26l's `convert_to_layered`. The dimensions panel shows "Unlock shapes (edit in the 3D modeller)" for a
  3D tensor deck (the preview) and an extrude3d deck (the modeller); a refusal shows in red under it and nothing is
  written. Once unlocked the Geometry tab shows the 3D modeller, and its Undo takes the deck back to the tensor deck
  and the preview (`SolidPage::record_undo`; `MainWindow::solid_edited` switches the page back).
- **Measured.** DC8 (hidden `[.unlock_library]`, Release): the planar MOSFET (3D) converts (s = 1.2 nm, max_size
  100 nm, estimate 119k nodes) and builds: 170,711 nodes, 285 s, its regions equal to the tensor deck's at 200 seeded
  points (points within 2 s of a box face left out). The bulk FinFET is refused by the 3D generator: "a patch edge
  meets an interface other than at a single trace in its face (unsupported)". The SOI FinFET is refused: its layout
  needs s = 0.025 nm. DC11: 12 of 12 mutants killed (one, the linspace spacing, after a check on a coarse x axis). DC12: `app/unlock` has no Win32 header. Targeted suites pass (`[unlock]`, the
  dimensions, geometry, modeller, solid, sweep and New device window files and model tests).
- **Rulings.** The automatic s comes from the rules alone (`largest_resolution`), not trial meshes
  (`largest_meshable_resolution` meshed the MOSFET for over 20 minutes); the build then meshes it. No node cap: the
  estimate is in the note. For an extrude3d deck Unlock does what the inspector's Convert does (both kept).
- **Limitations (recorded).** Neither library FinFET unlocks (the generator's rules, Unit 64). An unlocked deck's
  mesh is far larger than its tensor mesh (MOSFET 3D: 170,711 nodes against 8,525). After Undo back to the preview,
  Redo is gone (the 2D page's history restarts there).

**Slice 4 record (2026-11-23, branch `app/device-templates`, commits after slice 2's, local only):** DT8, on the
owner's word ("after slice 2 is done, start slice 4 (choose recommended for decisions)").

- **Built.** `solid::outside_domain` (the first solid past the domain: region, face and distance), `clamp_to`,
  `clamp_offset`, `clamp_scalar` (bisection); `geo::outside_domain`, `geo::domain_holds`. 3D modeller: points on the
  work plane are held on the domain (a box's base, a cylinder's centre, a prism's vertices), heights and radii are
  held, and so are the move gizmo's offsets (axis, plane, move, an array's farthest copy). Face drags, Extrude and
  Inset depths stop at the last offset inside, found by bisection, and vertex drags are held. Any other drag frame
  that would leave a solid outside is not taken ("Held at the domain: ..."). `commit_model` refuses an edit that
  takes a solid out (a typed size, a smaller domain) with the reason ("Refused: region 'ox' reaches 0.3um past the
  domain's z_max face (a solid stays inside the domain).") and leaves the deck unchanged; a deck already past its
  domain is not held. The inspector's position and size fields show and take references (`$W`, `$W + 1um`); a number
  typed over a reference replaces it. "Name this value…" on the fields' context menu (a box or the domain) adds
  `L<k> = <value>` and writes `$L<k>` in the field, one undo. 2D editor: snapped points are held on a domain that is
  given (a tensor deck's, a layered deck's `domain` line); a commit taking a shape out is refused the same way.
- **Measured.** DC7 (`[guardrails]`, Release, WARP): 500 seeded random drags of Box, Prism, Cylinder and the move
  gizmo, 383 of them written, none leaving a solid outside (in the model or the deck), 14.4 s; a typed size past the
  domain is refused with the reason and the deck byte-identical. DC1: one 26m check changed: MC6's Inset drags 0.1 µm
  up from a top at 0.55 µm in a 0.6 µm domain, now held at 0.05 µm (the check said "depth 0.1"; it now says
  "depth 0.05"); recorded, not loosened. DC11: 10 of 13 mutants killed (two after a test of a box drawn and moved past the faces; three rewritten after their first form did not compile, two of them then killed). Survivors: the cylinder radius clamp (no test draws a cylinder at a face), the move/plane offset clamp (the gizmo tests drag an axis arrow, clamped on its own path) and the frame guard (behind the clamps). DC12: the new model code has no Win32 header.
- **Rulings.** What is typed in a field is what is written: a number over `$Lg` replaces the reference (MT8's ties
  still follow drags). Name this value takes the next free `L<k>`, as N does while dragging. In 2D a drag frame is not
  held mid-drag; the release is refused if a shape would leave the domain. A generated 2D deck without a `domain`
  line holds nothing (its domain is its shapes' bounds).
- **Limitations (recorded).** Contact patches are not clamped while moving (they stay on their face). Name this value
  is for boxes and the domain, not prisms. The menu replaces the edit's own context menu (Cut, Copy, Paste kept; Undo
  and Select all not).
- **Full suites (owner "run Full Release/Debug suites", 2026-11-23), at `ed29b4d`:** Release 301/301 pass (4038 s; one
  thread per test, `-j 6`; `verify_full.tvb04` 1413 s alone on one thread, within its earlier 1284-1880 s). Debug: the
  build first failed on two C4702s (a `FAIL` before a helper's return in the new_device and sweep window tests; now
  `FAIL_CHECK`). The suite then ran to 251/257 before being stopped on the owner's word: `app.solid` failed on MC9's
  1 ms tied drag frame (15.1 ms in Debug; now checked in Release only, as KG6 and SE5, the limit unchanged; Release
  0.668 ms), `app_window.dimensions_window` failed on a race (the open's build still running in Debug at the test's
  edit failed on the broken deck and went to the Deck tab; the test now waits for it), and `app_window.unlock_window`
  ran over 85 minutes (the extrude3d case unlocked the library FinFET, whose 0.5 nm layers had the live check mesh the
  converted deck; it now uses a small extrude3d deck, 3.7 s in Debug, the FinFET's conversion kept in `unlock_test`).
  With the fixes (b5163ad) the six tests left (those three, `app.run_controller`, `app_window.self_test` and
  `app_window.analysis_window`) pass in Debug, and the changed tests pass in Release. Limitation: unlocking the library FinFET in
  the app sets off that long background mesh build (the window stays responsive).

### Unit 64 — interface pieces the 3D generator refuses: holes, domain edges, slanted walls (proposed 2026-11-21; on `mesh/interface-pieces`)

**Owner request (2026-11-21):** asked what limits make the 3D generator refuse drawn fillets, cylinders and arrays
(Unit 26m's biggest risk) and how to fix them; "yes" to proposing this unit, holes first. The branch starts from
`app/solid-modeller` `a8ca49f` (Unit 26l, PR #111 open), whose layered 3D decks are the fixtures.

**Where things stand.** The polyhedral generator (Unit 53) meshes each interface piece by mirrored node pairs at
± s/2 (s = `interface_resolution`), so the dual faces lie in the interface and the P4 certificate (Unit 16) holds.
Its layout step (`src/mesh/complex_interfaces.cpp`) refuses:
- **holes:** a piece's boundary is walked as one loop; a second loop is refused ("an interface piece has holes",
  line 326). `pair_set` takes that one inset polygon. A nanowire through a gate, or a plug through a layer, gives
  such a piece;
- **domain edges and corners:** a piece edge whose two boundary facets lie on different domain faces is refused
  (line 373), and so is a piece reaching a domain corner (line 524). Fins, fillets and arrays flush with the domain's
  sides give these;
- **a plane crossing a slanted wall** (26l SD6, recorded): the 4 s clearance check (line 617) skips touching
  non-parallel pieces and exactly coplanar ones. A slanted wall split by a plane is two pieces folded by about an
  ulp (26l SQ4): they touch, are not exactly coplanar, and are parallel, so their distance 0 is refused. This
  reading is a hypothesis from the code, confirmed or not in slice 0.

**Standing conditions:** every deck and complex that meshes today meshes to the same bits (the new paths only run
where the old code refused). The P4 certificate is required and unchanged. The 60-degree sector rule (U12, the
owner's ruling of 2026-10-14) and the 4 s rule between unrelated interfaces are kept. Criteria fixed before code.
Targeted Release tests; full suites and Debug only on the owner's word; `-j 6`, 6 threads. No new dependency.
Nothing pushed or merged without authorization.

**Left out:** sectors under 60° (needs the measurement U12 asks for, and the owner's word to loosen a frozen rule);
oblique domains (layered domains are boxes); pinched pieces (two loops meeting at a vertex).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| IP1 | 0 | **Reproduce and pin.** A fixture per refusal, from layered 3D decks: a square and a 16-sided nanowire through a gate (hole), a fin flush with the domain's side and a fillet running into it (domain edge), a fin at a domain corner, a slanted fin crossed by the STI top (fold). Each refused today with the message named above; the fold hypothesis checked and the record corrected if wrong. | Tests only. |
| IP2 | 1 | **Holed pieces.** A piece's boundary is any number of loops: the outer one and its holes (orientation from the piece's normal). The inset is applied to every loop (holes grow by 1.5 s); the graded feature pairs run along every loop's edges; `pair_set` triangulates the inset domain with its holes. Hole edges are feature lines like any other (junctions with the piercing solid's side faces, sectors checked as now). | IQ1. |
| IP3 | 2 | **Domain edges and corners, right angles only.** A trace that reaches a domain edge where two box faces meet at 90° is straddled into both faces: the pairs are mirrored across each face plane, so by the same mirror symmetry as U12's exact right-angle bend the dual faces stay planar there. A domain corner (three faces at 90°) likewise. Any other edge angle stays refused. | Exactness argued as U12's 90° case; checked by IC3. |
| IP4 | 3 | **Folded coplanar pieces.** Two pieces that share an edge, have the same two regions on the same sides, and whose normals agree within 1e-12 are one 180° junction (skipped by the clearance check, as exactly coplanar pieces are now); the feature between them gets no corner nodes. | Only for 26l's ulp folds. |
| IP5 | — | **Commits:** this proposal, a commit per slice, ARCHITECTURE (the mesh chapter's generator rules, 11.3) and the records. The owner may stop after any slice. | |

**Criteria (fixed before code; a target missed is recorded as a limitation, never loosened):**

| # | Criterion |
|---|---|
| IC1 | **Nothing else moves:** every example and test deck, and every complex in the mesh tests, meshes to the same bits as before. |
| IC2 | **Certified:** every IP1 fixture now meshes, passes every import check, and carries a valid P4 certificate. |
| IC3 | **Accurate:** each fixed fixture's equilibrium (potential and carrier densities at the nodes, interpolated to 200 sample points) agrees to 1e-6 relative (1e-9 V absolute where the potential is near zero) with a reference: the same device meshed by the swept generator (Unit 54) where it is an extrusion (the nanowire, the flush fin), else the same layout moved 10 s off the domain edge (a domain-edge case) or cut exactly coplanar (the fold). |
| IC4 | **Still refused:** sectors under 60°, unrelated interfaces under 4 s, pinched pieces, non-right domain edges: the same messages as before (table). |
| IC5 | **Random holes:** 200 seeded layouts with 1 to 4 holes per piece (boxes and prisms piercing slabs) mesh and certify, or are refused only by a rule IC4 keeps (counted). |
| IC6 | **Cost:** a holed fixture's node count ≤ 1.2 × that of the same layout with the hole's solid made of the slab's region plus the hole's own interface (recorded, not a gate, if missed); generator time recorded. |
| IC7 | **Mutants:** at least 8 per slice caught (e.g. a hole loop not inset; a hole's orientation reversed; pairs mirrored across one face only at a domain edge; a non-right edge accepted; the fold tolerance applied to non-touching pieces). |

**Questions (each with a recommendation; ruled by the standing "choose recommended" rule, 2026-11-11):**

- **IQ1 — holes.** (a) Native loops: the inset, the graded pairs and `pair_set`'s triangulation take holes. (b) Seams:
  cut each holed piece into simple pieces along added edges, meshed as now (the seams become 180° junctions with
  extra pairs). **Recommended (a)**: no artificial feature lines, fewer nodes; (b) only if (a) cannot certify.
  Ruled (a).
- **IQ2 — order.** (a) Holes, then domain edges, then folds (by how often a drawn device hits them: every
  gate-all-around device has a hole). (b) Folds first (smallest). **Recommended (a)**. Ruled (a).
- **IQ3 — the reference in IC3.** (a) The swept generator where the device is an extrusion, else a shifted or
  exactly coplanar layout. (b) A tensor grid. **Recommended (a)**: same physics path, different mesher. Ruled (a).

**Risks:** the domain-edge argument assumes mirror symmetry, which grading toward a corner may break (U12's
measurements found bisector nodes made a right-angle corner worse), so IC2 and IC3 are the gates; `pair_set` with
holes is new triangulation code, so IC5's random layouts guard it.

**Slice 0 results (2026-11-21; IP1; `be45c2d`):** tests only, `tests/mesh/interface_pieces_test.cpp` (ctest
`mesh.interface_pieces`, 4 cases, 25 assertions, 7.2 s Release); the generator unchanged. Fixtures are layered 3D
layouts (`layered_partition3d`) in a 1 x 2 x 1 µm oxide domain, s = 0.02 µm unless named.
- **Holes: confirmed.** A square wire (0.2 µm) and a 16-sided wire (radius 0.2 µm, s = 0.01 µm so its traces clear
  6 s) through a gate: the gate's faces are pieces with the wire as a hole, refused "an interface piece has holes".
- **Domain edges and corners: the record was wrong.** A fin flush with the domain's x = 0 face, the same fin under a
  gate flush with the top and both x faces, and a block in the domain's corner touching three faces all **mesh
  today**. A filleted flush fin (radius 0.2 µm, four segments of 0.078 µm) meshes at s = 0.01 µm and is refused at
  0.02 µm by the 6 s trace rule (a sizing rule, kept). The reason: an axis-aligned piece in a box domain never has
  an edge on a domain edge or a vertex at a domain corner (its plane would have to hold the edge); where it reaches
  one it crosses it, which Unit 53 handles. The refusal "an interface ends on a domain edge" is reached only by a
  slanted wall ending on the edge (pinned, a wall from (z, x) = (1, 0) to (0.4, 0.2)); with the faces at 90° such a
  wall leaves a wedge of at most 45° on one side, so the kept 60-degree rule (U12) refuses it if this rule did not.
  A domain corner is reached only with an edge (26l's slanted walls are all parallel to one axis), so "a domain
  corner" is not reached on its own. **IP3 (slice 2) therefore has no target while U12 stands: dropped**
  (recommended; ruled by the standing rule of 2026-11-11; dropping a slice loosens nothing). Reviving it needs the
  owner's word on U12 for sectors under 60°.
- **Folds: mechanism confirmed, scope corrected.** 26l SD6's slanted FinFET is refused "two interfaces within 4 s".
  Its touching, not exactly coplanar wall facets are parallel to 1e-15 (the rule's tolerance 1e-12): 8 pairs of one
  region pair (a wall split into pieces of its own pair, a 180° bend) and 6 of different pairs (fin | oxide against
  fin | gate where the gate's bottom crosses, a 180° junction). IP4's "the same two regions on the same sides" is
  too narrow: slice 3 takes both kinds.
- **IC7:** 10 mutants of the rules the fixtures reach (hole check off; parallel tolerance 0 and 1e-17; touching
  parallel pieces skipped; coplanar neighbours not skipped; domain-edge rule off and inverted; 6 s rule at 3 s and
  off; the 4 s rule at 0): all 10 caught.
- **IC1:** the generator is unchanged; `mesh.generate_complex_interfaces`, `mesh.layered3d` and
  `mesh.layered3d_slanted` pass. Full suites and Debug not run (owner's word needed).

**Slice 1 results (2026-11-21; IP2, IQ1 a):** holed interface pieces meshed natively, in
`src/mesh/complex_interfaces.{hpp,cpp}` only.
- **Code:**
  - **Loops:** a piece's boundary is walked into every loop; the outer one is the loop of largest area in the frame
    (it fixes the frame). The rest are holes, clockwise, in `hole_loops`, `holes`, `local_holes` and `inset_holes`.
    The first loop still starts at the lowest vertex, so a piece without holes takes the old path unchanged.
  - **Checks on hole edges:** sectors (U12), feature-line length (6 s) and separation (4 s). The piece's distance
    (closed polygon less the holes' interiors) and its clearances run over the hole rings too. A hole meeting the
    boundary would be refused as internal (it cannot: it would pinch the outer loop).
  - **Inset:** each hole grows by 1.5 s, with the same offset to the left of its edges. It must stay simple,
    clockwise, keep its edges' directions, lie inside the inset polygon and stay clear of the other insets, else
    "too narrow for its 1.5 s inset".
  - **Pair set:** a holed inset domain is cut along lines v = c into simple polygons of one region (`cut_at_holes`).
    There is one line per hole, at the middle of the widest gap between the rings' vertex v values inside the hole's
    v range, so no vertex lies on a cut. The cut parts inside the domain (by parity) become edges. The faces are
    traced left of the rings and of both directions of each cut: arriving on a ring edge at a cut, turn onto it;
    arriving on a cut, continue on the ring. Every hole is cut, so every face is simple. The 2D generator meshes the
    tiling, and the hole edges' graded lattice is placed as on the outer edges. The cuts are not 3D features.
- **IC1 met:** five layouts without holed pieces (the flush fin, the flush fin under a gate, the corner block, the
  filleted flush fin at s = 0.01 µm, the uncrossed slanted fin) mesh to the same bits as at `be45c2d`. The pinned
  FNV-1a digests over points, cells and node regions were taken from the old generator. One existing test moved,
  as intended: Unit 53's FC2 section "an interface piece with a hole" (`generate_complex_interfaces_test.cpp`)
  pinned the refusal, and now requires that extruded slab with a hole to mesh and certify. The public header's
  interface rules (`generate.hpp`) no longer list holes as refused.
- **IC2 met:** the square wire (7,700 nodes) and the 16-sided wire at s = 0.01 µm (23,151) mesh, pass the import
  checks and carry a valid P4 certificate (`certifies`, flat deviation <= 1, corner deviation <= s). Their volumes
  sum to the domain's to 1e-12, and their digests are pinned.
- **IC3 missed (recorded, thresholds unchanged; the hidden case `[.interface_pieces_accuracy]` is `[!shouldfail]`,
  as Unit 53's FC8).**
  - **Set-up:** a p-type (1e17) silicon wire through an HfO2 block in oxide, with an ohmic body on the wire's y = 0
    end and a metal gate (4.4 eV, 0.5 V) on the four faces along the wire; equilibrium Poisson; 200 seeded points.
  - **Reference:** the swept generator's mesh of the same device, from two rings of trapezoids cut on the
    diagonals. A grid of cells failed the 2D generator, first on 6 s, then on P4 where internal edges end on an
    interface. With the gate on the top face only, gate and oxide cells became separate 2D regions; with the gate
    over the ends too, the gate's edge met the body's at the wire.
  - **Measured:** holed against swept, potential max 64 mV, rms 12.8 mV (of 0.71 V). The swept mesh against itself
    at h/2 and s/2 differs by 43 mV max, 13 mV rms, so 1e-6 is far below what two meshes agree to.
  - **Finding:** the two generators converge to different potentials. At (0.3, 1.0, 0.5) µm, in the HfO2, the
    holed mesh gives -0.174 / -0.165 V at h, h/2, and the swept -0.206 / -0.233 / -0.246 V at h, h/2, h/4. The gap
    remains with the HfO2 at the oxide's permittivity, and with no block at all (no holes: the polyhedral path that
    slice 1 leaves bit-identical). So it is not slice 1's; the contacts select the same faces on both meshes, and
    both meshes' numbering agrees. Which mesher, or this swept set-up, is off is an open question, outside this
    unit.
- **IC4 met:** two wires 2.5 s apart, a wire 2.5 s from the gate's side ("two unrelated interface features within
  4 s"), two wires touching along an edge ("pinched vertex") and a diamond wire with 30-degree corners (U12) are
  refused with the same messages as before.
- **IC5 met:** 200 seeded slabs pierced by 1 to 4 boxes and prisms (4 to 6 sides), 464 holes in all: 200 mesh and
  certify, none refused; 468 s in Release. The first 20 run in the visible case; all 200 in the hidden
  `[.interface_pieces_random]`.
- **IC6 met (recorded):** the square wire takes 7,700 nodes against a reference of 7,515 (the gate alone, plus what
  the wire alone adds to the empty domain): 1.025, under 1.2. The generator took 1.9 s.
- **IC7 met:** 12 mutants; 9 caught. These were: hole not inset (0 s); hole orientation reversed; distance
  ignoring holes (caught once the holed digests were pinned); cut parity shifted; the cut never taken (nested rings);
  holes not passed to the pair set; the cut at a vertex; the smallest loop taken as outer; only the first hole cut.
  3 are equivalent: hole loops left out of the side checks and of the sector check (every hole edge is also an edge
  of the piercing solid's pieces, whose loops carry it), and hole rings left out of the clearances (`distance()`
  already measures to the holes).

**Slice 3 results (2026-11-21; IP4, both kinds as slice 0 found):** folded slanted walls, in
`src/mesh/complex_interfaces.cpp`, and (owner's ruling, below) the P4 bound for oblique non-convex pieces in
`src/mesh/interfaces.cpp`.
- **Code, the generator:**
  - **Fold:** two interface facets on one edge, not exactly coplanar, with normals parallel to 1e-12 (the clearance
    rule's own tolerance), and lying on opposite sides of the edge (a 180° bend, not a 0° one).
  - **Same pair:** folded facets join one piece, as coplanar ones do.
  - **Different pairs:** pieces sharing a fold edge are skipped by the 4 s clearance check, as exactly coplanar
    pieces are (a 180° junction). Touching parallel pieces without a fold edge are still checked.
  - **Corners:** a loop vertex on a fold edge whose two edges are straight to 1e-12 (sine) is not a corner (corners,
    traces, sides). Elsewhere collinearity stays exact.
  - Every layout with a fold was refused before ("within 4 s" at distance 0), so the new paths run only there.
- **Finding: a second refusal behind the fold.** With folds taken, SD6's FinFET was refused by P4b ("a corner
  zone deviates beyond the resolution"). The exactly coplanar version on a dyadic grid, which has no fold, was
  refused the same way. So it is not the fold.
  - Probes: the refusal needs a fin | oxide wall piece wrapped around the gate (L or U shaped, so non-convex) on the
    slanted wall. A gate along the whole fin, no gate, an uncrossed gate, or an upright fin under the same gate all
    certify.
  - The bound was 1.08 to 1.35 s, against s. The mirrored pairs put the faces in the plane to rounding, so the
    excess is the certificate's own: for oblique non-convex pieces it used only the decomposition's min over
    triangles of max over V, too large when V straddles a diagonal (its comment says so).
- **Owner's ruling (2026-11-21), asked because the standing conditions say the P4 certificate stays unchanged:**
  tighten the bound, acceptance limits unchanged.
  - **New bound:** each face vertex x maps to its orthogonal foot on the plane through v[0], taken in the exact
    projection the ear clipping uses (dropping the normal's largest component n_drop):
    q(x) = project(x) - h project(n), with h = (x - v[0]) . n.
  - **Soundness:** q is affine, so the hull of q(V) is q of V's hull. If it lies in the projected polygon grown by
    an offset, every point of the face lies within max |h| plus the offset's reach / |n_drop| of the polygon (the
    projection shortens in-plane lengths by at most |n_drop|).
  - **Offsets:** 1e-14 X first (X the largest coordinate magnitude: the feet are rounded by a few ulps of X, well
    inside P4a's 1e-13 X), then 1e-6 s up to s / 4 as for axis-aligned pieces. The exact polygon is not tried
    (q is rounded).
  - **Kept bound:** the smaller of this and the decomposition's. A first try at 1e-6 s alone passed P4b and failed
    P4a: the offset, about 1e-12 cm, was above P4a's tolerance.
  - **Effect:** a smaller upper bound, never a larger one. Reported deviations of existing oblique non-convex
    pieces can only fall, and no test pins one (upper bounds only). Mesh bits do not move.
- **IC1 met:** the five pinned digests (slice 1's) unchanged. 26l's SQ4 case
  (`layered3d_slanted_test.cpp`) moved as intended: it pinned SD6's "within 4 s" refusal and now requires the mesh
  and certificate. The whole `nitcad_mesh_test` passes in Release (180 cases, without the 200-layout IC5 case,
  20 min). Full suites and Debug not run (owner's word needed).
- **IC2 met:** SD6's FinFET meshes and certifies (7,136 nodes, pinned). The dyadic FinFET certifies both exactly
  coplanar (k = 1, 9,007 nodes) and folded (k = 1 - 1e-9: 8 same-pair and 6 different-pair folds, 8,927 nodes),
  each pinned. The exact one had been refused by P4b before the bound changed.
- **IC3 missed (recorded, thresholds unchanged; the hidden case is `[!shouldfail]`):**
  - **Set-up:** p-type silicon (1e17), region 2 as HfO2, ohmic body on z = 0, metal gate (4.4 eV, 0.5 V) on z = 1;
    equilibrium Poisson; 200 seeded points scaled by k.
  - **Measured:** folded against exact cut, potential max 26 mV, rms 6.0 mV (of 0.71 V).
  - **The floor:** the meshes cannot be the same. The generator's node placement is not continuous in the
    geometry: the uncrossed FinFET, with no fold at any k, moves from 8,371 to 8,200 nodes under the same 1 - 1e-9
    scaling, and by 3 nodes under 1e-12. Its solutions differ by 23 mV max and 4.9 mV rms, close to the folded case.
  - **Reading:** the fold adds little over any ulp-level change of the layout, and 1e-6 against an exactly coplanar
    cut is not reachable by a generator that places nodes discontinuously. IP4's "no corner nodes" holds at the
    layout level (mutants 5, 13 and 14 below), not as equal meshes.
- **IC4 met:** a gate block 1.5 s above a slab (parallel, overlapping, not touching) is still refused "within 4 s";
  the slice 1 refusals are unchanged.
- **IC7 met:** 16 mutants, 12 caught.
  - Caught: the fold tolerance at 0; the opposite-side test inverted; the fold skip removed; the fold skip for same
    pairs only; the same-pair union off; fold vertices kept as corners (all of them, one end of each fold edge, or
    by a 1e-17 straightness tolerance: the last three abort on the contract `std::isfinite(x)` in
    `exact_integer.cpp`, recorded below); the oblique bound off; its first offset at 1e-6 s; the largest offset
    tried first; the lift by the normal's smallest component.
  - Equivalent:
    - The fold skip widened to every parallel pair: overlapping parallel pieces under 4 s always leave a non-parallel,
      non-touching pair or a domain-edge crossing under 4 s too, refused by the kept rules.
    - Tolerant straightness at every vertex: no fixture has a near-straight vertex off a fold.
    - The reach not lifted by 1 / |n_drop|, and the feet without the h correction: at 1e-14 X and h at rounding
      neither changes an outcome, though both are needed for the bound to be rigorous.
- **Recorded, outside the unit:** a corner at an ulp angle (the mutants above) aborts the generator on a contract
  instead of being refused with a message. Slice 3 never forms one; any rule that does would need a guard.
- **Runs:** targeted Release, `-j 6`, test shards at most 6 at once.
- **Full Release suite (owner's word, 2026-11-21; at `d79bdb9`, `ctest --preset release -j 6`):** 278 of 279
  passed in 825 s (1 skipped: the PARDISO MOSFET case, as always). The one failure was slice 1's, unseen then (no
  full suite): `app_window.solid_check_window` used a box punching a hole in the oxide's interface plane as an
  invalid deck. The hole meshes now, so the deck is refused by the next rule ("an interface feature line is shorter
  than 6 s", the box's 0.1 µm sides) with the marker on the box's side (z = 0.45 µm, not 0.4). Still one
  diagnostic naming the interface: the row's name, comment and expected z were updated, and it passes.
- **Full Debug suite (owner's word, 2026-11-21; at `07372d2`, `ctest --preset debug -j 6`, a fresh Debug build,
  /W4 /WX clean):** 242 of 242 passed in 1,048 s.

### Unit 62 — density-gradient quantum confinement, opt-in (proposed 2026-11-19; on `physics/density-gradient`)

**Owner request (2026-11-19):** "Plan for Density-gradient confinement (Commercial TCAD feature type)". This is the
open half of ARCHITECTURE.md 11.3 "Quantum confinement": Unit 47 built MLDA, and density gradient was left as
"planned, no unit named". Docs only: no code until the questions below are ruled and the owner says start.

**Why density gradient and not more MLDA.** MLDA is a local factor on the densities of states. It needs the distance
to an interface and has no current-carrying physics of its own. Density gradient (Ancona and Iafrate 1989; Ancona and
Tiersten 1987; on unstructured grids Wettstein, Schenk and Fichtner, IEEE TED 48 (2001) 279) is the form commercial
simulators ship as their default confinement model. It adds a quantum potential Λ per carrier to the band edge
seen by the transport, solved with ψ, n and p. It confines at any interface shape with no distance search: fin
corners, the wrapped gate stack and the curved Unit 53 geometry included. It also acts along the channel: the source
barrier and short-channel quantum effects on the current. MLDA cannot do this.

**Standing conditions.** Opt-in. With the model off, everything is bit-identical to `main`: residuals, Jacobian
patterns and values, run identities, result files, and the Unit 40/59 traces and decks. Criteria are fixed before
code and never loosened or tuned after a run (TVM2). The coefficient is not fitted to pass a gate. Release and Debug
suites run with `/W4 /WX`, and only on the owner's word. Commits stay local, with no push, PR or merge without the
owner's word.

**Model (scaled by k T and V_T as 6.1).**
- Densities: n = n_ie exp(ψ − φ_n − Λ̂_n) and p = n_ie exp(φ_p − ψ − Λ̂_p) (Boltzmann; Fermi–Dirac in QP7), with
  Λ̂ = Λ / k T. Λ_n > 0 raises the electrons' effective conduction band, and Λ_p > 0 lowers the holes' valence band.
- Quantum potentials: Λ_n = −(γ_n ħ² / (6 m_n)) ∇²√n / √n, and the same for Λ_p with p and m_p. Since
  ∇²√n / √n = ½ ∇² ln n + ¼ (∇ ln n)², this equals −(γ_n ħ² / (12 m_n)) [∇² ln n + ½ (∇ ln n)²]. With
  β_c = γ_c ħ² / (6 m_c k T) (a length²), Λ̂_c = −β_c ∇²√c / √c. γ_c is a dimensionless prefactor whose
  theoretical value (and the mass it pairs with) is **not yet verified against the primary sources** (QP0).
- Transport: the Scharfetter–Gummel fluxes take ψ − Λ̂_n (electrons) and ψ + Λ̂_p (holes) as the potential across
  each edge. Poisson's charge, mobility and generation driven by field or current read the densities, which already
  carry the confinement.
- Recombination and traps (QP13) are **not** density-only. At DG equilibrium n p = n_ie² exp(−Λ̂_n − Λ̂_p), so terms
  that compare n p with n_ie², or use n₁ and p₁, must see the effective band edges, or they generate at zero bias.

**Design (QP0–QP15).**

| ID | Decision | Reason | Limitation |
|---|---|---|---|
| QP0 | **Source check, before any code (the first step of slice 0).** Read and cite, with equation and page: Ancona and Iafrate, PRB 39 (1989) 9536; Ancona and Tiersten, PRB 35 (1987) 7959; Ancona, "Equations of state for silicon inversion layers", IEEE TED 47 (2000) 1449; Wettstein, Schenk and Fichtner, IEEE TED 48 (2001) 279. Record: (i) the coefficient's exact form and its theoretical value (Ancona's r factor; the high-temperature value against the pure-state value), which fixes the default γ of QP6; (ii) the mass each source pairs with it for (100) silicon electrons and holes, which fixes QP6's default masses; (iii) the oxide boundary condition each uses, which confirms or revises QP5; (iv) the published DG-versus-Schrödinger–Poisson agreement for inversion layers (charge, centroid, threshold or C–V), which sets QV6's bounds. The bounds are then frozen in a DECISIONS commit before slice 1. If the sources cannot be obtained, slice 0 stops and the owner is asked. | The coefficient differs between sources by a constant factor, and the plan was written without the papers at hand (2026-11-19: the web copies tried were unreadable). TVM2 needs the default and the bounds fixed from a cited source before any run, not after. | Until QP0 is done, γ = 1 and QV6's numbers below are placeholders, not decisions. |
| QP1 | Switch: `[models] quantum_confinement = density_gradient` (the existing key, with values `off`, `mlda` and now `density_gradient`, so MLDA and DG exclude each other by construction). Keys: `dg_mass_n`, `dg_mass_p` [m0] and `dg_gamma_n`, `dg_gamma_p` (defaults from QP6). The run record names the model, the masses, γ and each insulator's penetration lengths (QP5). | One confinement model at a time. A deck that says `density_gradient` today is refused (`tests/deck/quantum_confinement_test.cpp`) and will now be accepted. | — |
| QP2 | Unknowns: Λ̂_n and Λ̂_p, built-in fields `quantum_potential_n` and `quantum_potential_p` (`FieldId`), interleaved after ψ, n, p (stride 5, as Unit 52). Held at 0 on insulator nodes (FD5, `FieldDomain::semiconductor`). Built-in table terms (terms.cpp and the equilibrium table), not the GG8 extension, which is steady-only. Energy balance also has stride 5, so every branch is chosen by the model flag (a new `SystemData::density_gradient_`, like `electrothermal_` after GG8) and the `DofMap`'s field ids, never by the stride. A mutant that tests the stride instead must fail (QV14). | Transient, small signal, equilibrium and results need the built-in path. GG8 shows the layout takes new fields without a second architecture. Unit 57 GG8 found `electrothermal()` defined as `m_ == 4`: the same trap at stride 5. | Per-iteration cost like energy balance (3.8–5.8× drift-diffusion, Unit 52); measured in QV13. |
| QP3 | Λ row at semiconductor node i (box method, divided by √c_i): V_i Λ̂_i + β Σ_j (A_ij / d_ij) (exp(½ (ln c_j − ln c_i)) − 1) − β B_i = 0, the sum over the node's semiconductor edges. B_i is the barrier term of QP5. Meshed oxide: the interface is at the midpoint of each semiconductor–insulator edge (mesh/interfaces.hpp; where MLDA puts it, Unit 47 QD2), not at a node. So each such edge ij contributes A_ij / (d_ij / 2 + x_p): the half-edge in the silicon in series with the decay length in the barrier, eliminating √c at the face. Lumped gate: the node is on the surface (d = 0, as QD3), so B_i = A_i^G / x_p with A_i^G the node's gate area. | The Wettstein et al. form: the ratio form never forms √c of a rounding-level density, and the row is linear in Λ. | DG in this form is the leading-order (ħ²) correction. Not a Schrödinger solution: QV6 measures the gap. |
| QP4 | Boundaries. Ohmic contacts: Λ = 0 (Dirichlet; the carriers are classical at a contact). Semiconductor–insulator interfaces: QP5. Every other boundary: natural (zero normal derivative of √c). | Standard DG boundary set. | Schottky contacts refused (QP10). |
| QP5 | Insulator interfaces: barrier penetration, analytic (Robin). In the barrier √c decays as exp(−x / x_p), with x_p = ħ / √(2 m_ox ΔE). ΔE is the band offset from the insulator band data of Unit 51 GD5 (affinities and gaps; electrons ΔE_c, holes ΔE_v), and m_ox the insulator's tunnelling mass. The band data come from the same lookup gate tunnelling uses: `band_data()` (gate_tunnels.cpp) takes the insulator's own data, or SiO₂'s tunnelling set for the built-in SiO₂ (which has none of its own, insulator.hpp). A lumped gate takes its stack's data, else SiO₂'s. That gives B_i of QP3. Refused only: a meshed insulator that is not the built-in SiO₂ and has no band data, or a stack whose permittivity is not SiO₂'s and has no band data (as gate tunnelling). | No oxide unknowns. It works for meshed and lumped gates alike. It is derived from the barrier, not fitted. Every existing SiO₂ deck is accepted unchanged. SiO₂ electrons (3.0 eV, 0.50 m0, GD5): x_p ≈ 0.16 nm, comparable with the 0.1 nm interface spacings in use, which is why QP3 keeps the half-edge. | Penetration solved in the oxide (Λ rows on insulator nodes) is the more exact form: a follow-up (QR3 b). |
| QP6 | Default γ and masses: those QP0 finds in the cited theory (the theoretical value, not a calibrated one), recorded with their source before slice 1. Placeholder until then: γ = 1 with MLDA's masses (0.98 and 0.49 m0, Unit 47 QD1). Deck-overridable (`dg_gamma_*`, `dg_mass_*`), and always in the run record. | No fitted default (TVM2). Commercial tools ship calibrated γ: here the owner may set one in the deck, and it is reported, never defaulted from a fit. | The theoretical γ may miss Schrödinger–Poisson by more than a calibrated one: QV6 records it honestly. DG's single isotropic mass is not the reference's per-valley m_l = 0.916 and m_t = 0.19 (QP14). That mismatch is part of what QV6 measures, and is stated in its report. |
| QP7 | Fermi–Dirac: the Λ row unchanged. Λ̂_n shifts the electrons' band edge in the density relation. In the scaled band diagram of `node_bands` (scaled_device.hpp: η = ψ + s, E_c = g_n − η, E_Fn = ln(n / n_ie) − ln γ_n(n) − η), DG makes it ln(n / n_ie) − ln γ_n(n) = η + E_Fn − Λ̂_n, i.e. n = N_c F_{1/2}(E_Fn − E_c − Λ̂_n). Holes: ln(p / n_ie) − ln γ_p(p) = −η − E_Fp − Λ̂_p. Under Boltzmann (ln γ = 0) this is the model line above. | The form commercial tools use under FD; it reuses Unit 14's degeneracy factor unchanged. | The ħ² coefficient is the non-degenerate one. Documented. |
| QP8 | Equilibrium: the equilibrium table gains the Λ rows, with carriers slaved by QP7's relation at E_Fn = E_Fp = 0 (Boltzmann: n = n_ie exp(η − Λ̂_n) and p = n_ie exp(−η − Λ̂_p), η = ψ + s; Fermi–Dirac through γ_n and γ_p as QP7) (unknowns ψ, Λ̂_n, Λ̂_p). The zero-bias state is the DG equilibrium. | Every sweep starts from a true equilibrium: currents at rounding (QV4). | Equilibrium becomes a 3-field Newton with DG on. |
| QP9 | Start and continuation: Λ̂ = 0 from the classical solution. If Newton fails at the first DG solve, a coefficient continuation runs. It is a new loop in solve, not Unit 28's step cutting, which steps bias, not β. It solves at β·k for k = 1/8, 1/4, 1/2, 1, each from the previous solution. A failed rung halves the gap to the last converged k, at most 4 times, then the point fails with the rungs listed. The rungs taken are in the run record. Later points continue as every sweep. Λ̂ is clipped ±max_update and measured like ψ: `System::update_size`, `apply_update` and Unit 59's `System::natural_size` are extended so Λ̂ enters as a potential (no asinh compression), tested by a case where only Λ̂ moves and each measure sees it. | Deterministic, and bounded at the first point only. | The ladder's use is reported in QV13, not hidden. |
| QP10 | Refused (invalid_input; the same message from the equilibrium and drift-diffusion assemblers): electrothermal, energy balance, thermionic interfaces, Schottky contacts, gate tunnelling, nonlocal BTBT, any semiconductor region not silicon (as QD6: the masses are global, not per region, so there is no exception), and the insulators without band data of QP5. Accepted: Lombardi, Masetti, Canali, local BTBT, impact ionization, SRH, Auger, interface traps and charges, both statistics, transient and small signal. | The refused ones read band edges or densities at interfaces in ways not designed for Λ (as QD8). Combining them is a follow-up, each with its own gate. | — |
| QP11 | Transient: the Λ rows carry no storage, so their mass-matrix rows are zero (index-1 DAE, as Poisson). Small signal: the Λ rows' partials enter the AC matrix. | Commercial parity: DG under every analysis. | — |
| QP12 | Results: Λ_n and Λ_p [eV] as nodal fields in `BiasPoint`, transient and small-signal results, result files, `nitcad plot` and the app's field list. The derived effective band edges E_c + Λ_n and E_v − Λ_p go in derived.hpp. The quasi-Fermi levels include Λ: in `node_bands` and every derived or plotted E_Fn, E_Fp, E_Fn = ln(n / n_ie) − ln γ_n − η + Λ̂_n and E_Fp = −η − ln(p / n_ie) + ln γ_p − Λ̂_p (QP7 solved for E_F), so they are flat at equilibrium (QV4). The classical bands E_c and E_v stay available as before. Unit 61 checkpoints carry Λ, and a restart with a different confinement model is refused. | The quantum potential is what users inspect in commercial tools. | Result files gain the fields only when the model is on. |
| QP13 | Recombination and traps see the effective band edges E_c + Λ_n and E_v − Λ_p. SRH (bulk and interface traps), Auger and radiative: n_ie² → n_ie² exp(−Λ̂_n − Λ̂_p), n₁ → n₁ exp(−Λ̂_n), p₁ → p₁ exp(−Λ̂_p), and the same substitution in Unit 14's Fermi–Dirac forms. Trap occupancy (charge and the trap's SRH) uses the same n₁ and p₁. A trap's declared level keeps its classical position (as QD9). The Λ partials of every such term enter the Jacobian (QV5). | Detailed balance: with these substitutions, every net rate is zero at DG equilibrium (n = n₁ e^{…} relations hold exactly), and the occupancy at equilibrium is the Fermi level's (QV4). Field- and current-driven generation (impact ionization, local BTBT) has no n p term and is unchanged. | Λ enters every recombination term's Jacobian: more entries per node. |
| QP14 | Verification reference: a 1D Schrödinger–Poisson solver in `verify/reference` (finite differences; a symmetric tridiagonal eigensolver of its own, Sturm bisection and inverse iteration as Unit 55 TVQ4 ruled, since `jacobi_eigen` takes n ≤ 16 and serves as its cross-check; silicon's 2- and 4-fold valleys with m_l = 0.916 and m_t = 0.19, a finite barrier with the insulator band data, Boltzmann and FD subband occupancy, self-consistent with Poisson). Electrons only: in a p-type MOS-C in inversion the holes are classical (bulk majority), so the reference tests Λ_n against physics, and Λ_p only through QV3–QV5 and QV8. A valence-band Schrödinger model (heavy, light and split-off holes) for an n-type accumulation case is a follow-up, not in this unit. Built through the Unit 55 §3.4 checklist: independence I2, class H. | DG has no closed form. SP is the physics DG approximates, and the reference commercial calibrations use. | 1D only (MOS-C). The FinFET is reported, not gated. Λ_p has no physical reference in this unit: stated in every QV6 report. |
| QP15 | Coarse meshes: across a depletion edge or junction on a mesh much coarser than λ, neighbouring densities can differ by decades, and exp(½ Δ ln c) makes Λ large and mesh-dependent (a known DG artefact on under-resolved meshes). No limiter or damping is added: that would change the model. QV15 measures it on the example decks, and the run record warns when \|Λ̂\| exceeds 10 (10 k T) at a node away from the interfaces: a semiconductor node with no semiconductor–insulator edge that is not a lumped-gate node. | Honest behaviour first. A limiter would be a model change needing its own ruling. | DG on the existing 5 nm example meshes may be unconverged in Λ: the warning says so, it does not hide it. |

**Criteria (QV1–QV15), fixed before code (QV6's bounds after QP0, before slice 1).**

| ID | Criterion |
|---|---|
| QV1 | Off: every existing test and pin passes unchanged. The Unit 59 NG1 trace set is byte-identical to `main`'s (only added cases new), the NG2 decks are 38/38, and run identities and Jacobian patterns are unchanged. |
| QV2 | β and x_p against long-double CODATA within 1e-14 relative; SiO₂ x_p printed. |
| QV3 | The Λ row alone, against analytic profiles: √c = exp(−x / 2L) gives Λ̂ = −β / 4L² exactly; √c = cosh(x / L) gives −β / L². Barrier term, both cases: (i) a meshed oxide, the interface at the half-edge, against the closed form of a 1D wall with the decay √c ∝ exp(−x / x_p) continued past the face; (ii) a lumped gate, the node on the surface. Both at spacings h = 0.4, 0.2, 0.1 and 0.05 nm, so h / 2 runs from above to below x_p, and the result converges as h → 0 with no jump between the cases. Error O(h²) on uniform and graded 1D meshes (observed order 1.9–2.1), and 2D/3D tensor meshes equal 1D within 1e-12. |
| QV4 | DG equilibrium (1D MOS-C in accumulation and inversion, the 2D MOSFET example, both statistics): every edge current within its rounding bound (as QC5); n and p against QP8's slaving relation (η = ψ + s; Fermi–Dirac through γ) within 1e-12 relative; every recombination and trap net rate (SRH, Auger, radiative, interface traps) within its rounding bound, with traps on; E_Fn and E_Fp (QP12) flat to rounding at every node. |
| QV5 | Jacobian against central differences (1D, 2D, 3D, Delaunay, lumped gate, meshed oxide; Boltzmann and FD; with Lombardi): largest relative column error ≤ 1e-6. Newton on the 2D MOSFET in strong inversion converges quadratically in its last three iterations. |
| QV6 | Against the QP14 Schrödinger–Poisson reference (1D MOS-C, 1.5 nm and 5 nm SiO₂, p-type 1e17 and 1e18, strong inversion at three gate biases, Boltzmann), with QP6's defaults. Three quantities: the inversion charge, the charge centroid, and the threshold shift (DG − classical) against SP's. **Bounds set by QP0** from the published agreement (Ancona 2000 and Wettstein et al. 2001): each bound is what those sources report for the theoretical coefficient, with its citation, frozen in a DECISIONS commit before slice 1. Placeholders until then, not to be used for a verdict: 10%, 0.3 nm and 25%. MLDA's values are reported beside them, and the per-valley mass mismatch (QP6) and Λ_p's lack of a reference (QP14) are stated. A miss is recorded as a fail, not tuned. |
| QV7 | Mesh convergence (QV6's 5 nm device): Q_inv at surface spacings 0.2, 0.1 and 0.05 nm. The last two differ by ≤ 1%, and by less than the first two. |
| QV8 | Limits: as x_p → 0 the interface density falls monotonically and the peak moves inward. As γ → 0 the solution tends to the classical one (difference ∝ γ, observed order 0.9–1.1). |
| QV9 | Symmetry: a y-uniform 2D MOS-C equals 1D within 1e-12; the 2D MOSFET example with x and y swapped gives the same I_d within 1e-9; a 3D extrusion equals 2D within 1e-12. |
| QV10 | Small-signal g_d and g_m against DC central differences within 1e-5; a transient long step equals DC within 1e-9. |
| QV11 | Every QP10 refusal, with its message from both assemblers; every QP10 acceptance runs. |
| QV12 | Results and restart: Λ fields in results files, `nitcad plot` and the app list; a checkpointed DG sweep restarts and finishes bit-identical to an uninterrupted one (as Unit 61 RG); a restart with the model changed is refused. The vocabulary pins grow with the new value and keys and are updated, not loosened: the deck's `[models]` list and the app's deck-completion list (`[models]` keys and `quantum_confinement` values; the pin that failed CI in Unit 61, 7bc4a15). |
| QV13 | Reported, not gated: the 2D/3D MOSFET examples' and the FinFET decks' V_th, swing and on-current with classical, MLDA and DG; Newton iterations; per-iteration cost against drift-diffusion; uses of the QP9 continuation. The FinFET decks run DG with gate tunnelling off (QP10 refuses the two together), and each report says so. |
| QV14 | Full Release and Debug suites (on the owner's word), zero warnings, hand mutants caught (a Λ sign, a missing barrier term, the barrier at the node instead of the half-edge, an unscaled β, Λ dropped from one SG flux, the equilibrium row left out, a branch on stride 5 instead of the model flag, the classical n_ie² left in SRH, Λ left out of E_Fn, Λ̂ left out of `natural_size`, and others). |
| QV15 | Coarse meshes (QP15): the 1D diode, 2D MOSFET and 3D MOSFET examples, unchanged meshes (5 nm surface spacing), with DG on. Gated: every point converges, or the failure is reported with the QP9 rungs. Reported: the largest \|Λ̂\| and where it is, the QP15 warnings, and the I_d change against classical. A refined-mesh rerun of the 2D MOSFET (surface spacing 0.1 nm) shows whether the coarse values were mesh artefacts. |

**Owner questions (QR1–QR8).**

| ID | Question | Options | Recommendation |
|---|---|---|---|
| QR1 | Unknowns | (a) Λ̂_n, Λ̂_p added (stride 5). (b) √n, √p replace n, p (every term changes). | **(a)**: existing terms stay; the Unit 52 precedent. |
| QR2 | Carriers | (a) Both, always. (b) Per carrier (`n`, `p`, `both`; stride 4 or 5). | **(a)**: one layout. Per carrier is a later option. |
| QR3 | Insulator interface | (a) Analytic barrier penetration (Robin, QP5). (b) Λ solved inside meshed oxides. (c) Hard wall. (d) Natural (no condition). | **(a)**: physical, no fit, works with lumped gates. (b) is a follow-up. |
| QR4 | Coefficient | (a) The theoretical γ and masses found by QP0 (placeholder γ = 1 and MLDA's masses), deck-overridable, reported. (b) Calibrated commercial-style defaults. | **(a)**: TVM2 forbids a fitted default. |
| QR5 | Equilibrium | (a) DG in the equilibrium table (QP8). (b) Classical equilibrium, DG from the first bias point. | **(a)**: a true equilibrium start. |
| QR6 | Combinations | (a) The QP10 refusals. (b) Also accept gate tunnelling now. | **(a)**: each combination gets its own gate later. |
| QR7 | Reference | (a) An in-house 1D Schrödinger–Poisson (QP14). (b) Published curves only. | **(a)**: an independent, reproducible reference (Unit 55 §3.4). |
| QR8 | Statistics | (a) Boltzmann and Fermi–Dirac. (b) Boltzmann only. | **(a)**: parity with MLDA and with commercial tools. |

**Rulings (2026-11-19):** QR1–QR8 ruled (a), as recommended, under the owner's standing instruction "always choose recommended" (2026-11-11). Not started: slice 0 waits for the owner's word.

**Slices.**
- **0. Sources and reference** (docs and verify only): first QP0. It records the coefficient, masses, boundary
  condition and QV6's bounds with citations, in a DECISIONS commit. Then the QP14 Schrödinger–Poisson solver and its L0
  tests (a square well's subbands against the analytic ones, the triangular well against Airy zeros, charge
  neutrality), plus the QV3 profiles. No product code.
- **1. Physics and deck:** parameters, β and x_p (QV2), deck keys, refusals (QV11), run record. The deck accepts
  `density_gradient` from this slice, so until slice 2 both assemblers refuse it (invalid_input, "density_gradient:
  not built yet"), tested: no deck runs with the model silently missing. Slice 2 removes the refusal.
- **2. Assemble:** fields, Λ rows, SG coupling, exact Jacobian, the equilibrium table; QV1, QV3–QV5, QV9.
- **3. Solve and results:** start and coefficient continuation, transient, small signal, results files, plot, app,
  checkpoints; QV10, QV12.
- **4. Validation and close:** QV6–QV8, QV13, QV15; QV14 on the owner's word; ARCHITECTURE 11.3 and 6.x as built.

**Risks.** (1) Cost: stride 5 is roughly energy balance's factorization cost; the FinFET runs may take hours (Unit 58
would help). (2) First-point convergence: QP9's continuation is the bounded fallback. (3) The theoretical γ may miss
QV6: it would be recorded, and the owner may then rule a reported γ. (4) The FinFET mesh-convergence check (11.3) is
still open, so FinFET numbers are reported only. (5) QP0 may not get the papers (the web copies tried on 2026-11-19
were unreadable): then slice 0 stops and the owner is asked. (6) Coarse meshes (QP15) may give large Λ on today's
example decks.

**Review (2026-11-19, owner "review the plan", then "fix them all").** Thirteen findings, all applied: (1) the barrier
term moved from the node to the half-edge for meshed oxides (QP3, QV3); (2) the band data from gate tunnelling's
`band_data()` lookup, so built-in SiO₂ decks are accepted (QP5, QP10); (3) non-silicon refused outright (QP10);
(4) the coefficient's source check made the first step (QP0, QP6, model); (5) QV6's bounds to be set from the cited
agreement before slice 1 (QV6); (6) the Λ_p reference gap stated (QP14); (7) coarse meshes (QP15, QV15); (8) the
mass mismatch stated (QP6); (9) the coefficient continuation specified (QP9); (10) the FinFET runs with tunnelling
off (QV13); (11) the vocabulary pins (QV12); (12) flags, not stride (QP2, QV14); (13) the FD relation written out
(QP7).

**Second review (2026-11-19, owner "now review again", then "yes").** Seven findings, all applied: (1) recombination
and traps see the effective band edges, so detailed balance holds at DG equilibrium (model, QP13, QV4); (2) the
quasi-Fermi levels include Λ (QP12, QV4); (3) the slaving relation in η = ψ + s with FD through γ (QP8, QV4);
(4) "away from the interfaces" defined (QP15); (5) `natural_size` extended and tested (QP9); (6) a "not built yet"
refusal between slices 1 and 2 (slices); (7) three more mutants (QV14).

**Slice 0 (2026-11-19, owner "start slice 0", then "proceed").**
- **QP0: blocked, owner supplying the papers.** None of the four primary sources could be read: Wettstein et al.
  2001 (IEEE, DOI 10.1109/16.902727, paywalled; abstract only), Ancona and Iafrate 1989 and Ancona and Tiersten 1987
  (APS, paywalled), Ancona 2000 (IEEE, paywalled); Ancona's 2011 review (DTIC copy HTTP 403) and the NASA and COMSOL
  PDFs fetched as compressed or scanned streams, not text. Secondary sources (a COMSOL tutorial's m0/3 DG mass) were
  not used: PC is not enough (Unit 55 §3.3). The owner chose to supply the PDFs and to build the reference meanwhile.
  γ, the masses and QV6's bounds stay placeholders; slice 1 does not start until QP0 is recorded.
- **QP14 built** (`include/NiTCAD/verify/reference/schrodinger_poisson.hpp`, `src/verify/reference/
  schrodinger_poisson.cpp`; class H, I2, provenance PD, derivation in the header): a symmetric tridiagonal
  eigensolver (Sturm count, bisection, inverse iteration with partial pivoting); the BenDaniel–Duke effective-mass
  equation by the box method with per-cell mass and potential; the self-consistent MOS capacitor (classical or
  quantum electrons, Boltzmann or Fermi–Dirac subband occupancy, classical holes), solved by outer iterations with
  a predicted density (its fixed point independent of the predictor), started from the classical solution.
- **L0 tests** (`tests/verify/reference/schrodinger_poisson_test.cpp`, 8 cases, 271 assertions, under 1 s Release;
  the whole `nitcad_verify_test` 178 cases pass): tridiag(−1, 2, −1) against its closed form (1e-14) and jacobi_eigen
  on n = 12; the infinite well against the discrete closed form and, at observed order 1.9–2.1, the continuum (also on
  a smoothly graded grid, 1.8–2.2); a finite well with a mass step against Brent on its transcendental equations; the
  triangular well against DLMF Table 9.9.1's Airy zeros (read 2026-11-19, PA); the classical Poisson part against
  `mos_capacitor.hpp` from accumulation to strong inversion; the quantum solution's Poisson residual, Gauss's law,
  orthonormality, and its changes with the grid halved (≤ 1e-3, ≤ 0.01 nm), 12 subbands instead of 8 (≤ 1e-4) and
  a 60 nm Schrödinger depth instead of 40 (≤ 1e-4).
- **Bounds changed after the first run (stated, not hidden).** Four L0 bounds set before the run were wrong a priori,
  while every observed-order check passed: (i) the discrete infinite well "1e-12 relative" ignored that eigenvalues
  are accurate to eps ‖A‖, not relative to E: now 4 eps ‖A‖ (≈ 4e-13 eV at 800 cells; observed 1.2e-13); (ii)–(iii)
  the finite and triangular wells' absolute bounds (1e-5, 1e-6 at the finest grid) underestimated the O(h²) constant,
  about (k h)² / 12 (observed 4e-5 and 1.7e-6): replaced by stricter checks, the Richardson value of the two finest
  grids against Brent within 1e-7 and against the Airy zeros within 1e-8; (iv) the classical MOS-C's single-grid
  bounds (1e-5 in Q_s, 1e-6 V) failed in accumulation at −1 V (1.4e-5, 1.2e-5 V): now observed order 1.8–2.2 between
  the grid and its halving, and the Richardson Q_s and band bending against `mos_capacitor.hpp` within 1e-7 (relative,
  and in V). None of these is a QV criterion; QV1–QV15 are untouched.
- **Reported** (1e17 p-type, 5 nm SiO₂, V_G − V_FB = 2 V, Boltzmann, (100) valleys 0.916/0.19 ×2 and 0.19/0.417 ×4,
  oxide 3.0 eV and 0.50 m0): N_inv 2.837e12 cm⁻² against 3.142e12 classically (−9.7%); centroid 2.167 nm against
  0.958 nm; 11 outer iterations; Fermi–Dirac N_inv 2.824e12.
- **Still in slice 0:** QP0 (on the papers), then the QV3 analytic Λ profiles (they need the DG row, so they move to
  slice 2, where the row exists).

**Ruling (2026-11-19, owner "choose recommended"):** slices 1–3 may start before QP0. The code does not depend on γ's
value: γ and the masses are inputs, their defaults marked as placeholders; QV6 gets no verdict until QP0 is recorded;
no criterion changes. Only the plan's order changes.

**Slice 1 (2026-11-19).** As built:
- `assemble::QuantumConfinement::density_gradient` (= 2) and `PhysicsModels::dg_mass_n`, `dg_mass_p`, `dg_gamma_n`,
  `dg_gamma_p` (defaults 0.98, 0.49, 1, 1: placeholders, `physics::dg_default_*`).
- `physics::dg_coefficient_cm2` (β) and `physics::barrier_decay_length_cm` (x_p) in physics/quantum_confinement.
  QV2 met: both against CODATA 2018 written out in the test within 1e-14 relative. SiO₂ electrons (3.0 eV,
  0.50 m0): x_p = 0.159373 nm; β_n(0.98 m0, 300 K) = 0.5013 nm².
- Deck: `quantum_confinement = density_gradient` and the four `dg_` keys; the bad-value message now names three
  values. App: the `[models]` vocabulary and sections list the value and keys (25 keys); the completion pin went
  from 20 to 24 (updated with the vocabulary, not loosened; QV12).
- Run record: with DG on, `models.quantum_confinement` (2) and the four `dg_` settings; MLDA's entries unchanged and
  written only with MLDA; DG's parameters are read only with DG on (identities tested).
- Refusals (QV11, `assemble::check_density_gradient`, called from `make_scaled_device`, so both assemblers give the
  same message): the masses and gammas, electrothermal, energy balance, nonlocal BTBT, gate tunnelling, Schottky
  contacts, a thermionic interface, non-silicon regions, a meshed insulator neither the built-in SiO₂ nor with band
  data, a lumped stack with no band data and a permittivity other than SiO₂'s. Everything accepted is then refused
  with "not built yet (Unit 62 slice 2)", tested on meshed and lumped SiO₂, a lumped stack with its own band data,
  SiO₂ with band data, and Lombardi + Fermi–Dirac + local BTBT.
- Run (Release, targeted, not the full suites): physics 5 cases, assemble 9, deck 4 (confinement filters), the
  whole `nitcad_app_test` (77 cases) and the completion window case: all pass.

**Slice 2 (2026-11-19, owner "start unit 62 slice 2 choose recommended").** As built:
- Fields `quantum_potential_n`, `quantum_potential_p` (kT, node, semiconductor). With the model on the unknowns are
  ψ, n, p, Λ̂_n, Λ̂_p (stride 5) and at equilibrium ψ, Λ̂_n, Λ̂_p (stride 3); every branch reads the slots, not the
  stride. Two table terms (src/assemble/density_gradient.cpp): the rows V Λ̂ − β B on the nodes, and the edge sums
  β g (√(c_j / c_i) − 1) (QP2), in both tables. On ohmic-contact and insulator nodes Λ̂ = 0 (QP4). At equilibrium the
  sums read the densities slaved to ψ and Λ̂ (QP8). The barrier terms B come from `make_scaled_device` (QP3, QP5):
  meshed oxide g d / (d/2 + x_p) on the semiconductor end, lumped gate A L_D^{2−D} / x_p, with x_p from the
  insulator's (or the stack's) band data. Gate tunnelling's band-data helpers moved to scaled_device so both read
  the same data. An insulator whose offsets to the semiconductor are not both positive is refused.
- `√(c_j / c_i) − 1` is formed through log1p and expm1 while the two densities are within a factor 1.5. The plain
  square root loses (a h)⁻² digits on a fine mesh, and measured 1.5e-12 between 2D and 1D at a h = 6e-4.
- SG coupling (QP6): electrons δ − ΔΛ̂_n, holes δ + ΔΛ̂_p, with `EdgeFlux::d_delta` for the Λ columns. Impact
  ionization gains the cross entries (n, Λ̂_p) and (p, Λ̂_n). The field mobility's driving field stays ψ's alone.
- Recombination (QP13): SRH (bulk and surface), Auger and radiative with E → E e^{−Λ̂_n−Λ̂_p}, n₁ → n₁ e^{−Λ̂_n},
  p₁ → p₁ e^{−Λ̂_p}, and every Λ partial (`physics::srh_recombination_shifted`). Interface traps take the
  semiconductor node's Λ̂ as constant over its half-cell: occupancy and trap SRH use the shifted n₁ and p₁, with Λ
  columns in the insulator's Poisson row.
- **Ruled (recommended): dopant levels keep their classical positions**, as trap levels do (QP13, QD9). Drift-
  diffusion ionization reads η_c + Λ̂_n and η_v + Λ̂_p, the densities' own band edges, with the Λ entries. At
  equilibrium the occupancy depends on η − log N_c, which the slaving leaves unchanged, so the equilibrium
  ionization is the classical one.
- Newton treats Λ̂ as a potential: in the update size, clamped to ±max_update, in the natural norm as a potential,
  and held on insulator and contact nodes. `System::state_from_quantum_potential(ψ, Λ̂_n, Λ̂_p)` gives the slaved
  state.
- **Ruled (recommended): the solvers refuse the model until slice 3.** The time step and the small-signal matrix
  are refused at the assemble layer by precondition. Every solver (bias sweeps, equilibrium, small signal,
  transient, traces, adaptation) returns invalid_input, "quantum_confinement density_gradient: the solvers are not
  built yet (Unit 62 slice 3)" (`assemble::density_gradient_solvers_unbuilt`). Both assemblers accept the
  combinations slice 1 refused as not built (strides 5 and 3 tested).
- Fixed on the way: `locate_surface` read energy balance's four rows from `m_ == 5`, which density gradient shares.
  It now excludes the model. energy_balance_ is set only after it runs, so the flag could not serve.

Results (Release, targeted, `nitcad_assemble_test` whole):
- QV3 (met): exp(−x/2L) and cosh(x/L) on the [0, 2 µm] bar, uniform and graded (x = s + s², ratio 3) at 20–160
  cells. Observed order 2.001–2.009 uniform and 1.950–1.991 graded. 2D and 3D extrusions equal 1D within 1.7e-14.
  Barriers: a profile √c = sinh((x + x₁)/L) that meets the wall's condition (Λ̂ = −β/L² exactly), at h = 0.4, 0.2,
  0.1 and 0.05 nm. Wall-node error (n): meshed −0.027, −0.036, −0.029, −0.019; lumped 0.86, 0.42, 0.21, 0.105
  (first order, h / 4x_p). Both decrease from 0.2 nm, once h/2 < x_p, and the first-order extrapolation reaches the
  closed form within 0.02 for both: no jump. Interior silicon: order 2 in both cases.
- QV4 (met in the assemble layer): DG equilibrium on the meshed and lumped 1D MOS-C (−1.5 V, 1 V) and the 2D MOS
  with its n+ pocket (1 V), both statistics, traps on:
  - the drift-diffusion residual on the slaved state is at most 4.3e-14 of each row's largest term (fluxes, net
    rates, trap terms and Λ rows);
  - the Boltzmann slaving n = n_classical(ψ) e^{−Λ̂_n} holds within 1e-12;
  - Newton's first update from it is at most 3.1e-12.

  The 2D MOSFET example, the terminal currents against their resolution and the flat E_F (QP12, results) need the
  solvers: slice 3.
- QV5 (met): house gate (largest relative column error, central differences), DD and equilibrium tables. 1D, 2D
  and 3D meshed oxide with traps, both statistics: ≤ 1.4e-7. Lombardi, impact ionization with field mobility,
  incomplete ionization and local BTBT: ≤ 1.4e-7. Lumped 1D and 2D, Delaunay, and a meshed oxide with band data:
  ≤ 1.1e-7.

  Two probe choices, measured:
  - The 2D and 3D MOS run without the n+ pocket. Across its abrupt junction the minority density jumps by 1e16, so
    the Λ rows there carry terms of 1e8 against an entry V. Their central difference is then rounding (up to 6e-3
    of the column), while the analytic entry matches the difference to 7 digits as the step grows to 1e-3.
  - The probes keep the densities smooth: no step-raised interface densities. The equilibrium probe is the
    classical solution with 1e-3 ψ noise, since 0.02–0.1 left Poisson rows of |f| 20–130 whose Λ columns were
    rounding at 1.0–1.8e-6.

  Newton on the 2D MOS (with the pocket) at 1 V from the classical equilibrium at 0.8 V, Λ̂ = 0: 23 iterations,
  last updates 6.4e-3, 2.9e-5, 5.9e-10 (order 2.0), then 2.8e-12 at rounding (6 iterations without the model).
- QV9 (assemble-layer part met): the lumped MOS-C at ±1.5 V, both statistics. 2D equals 1D within 6.2e-15, 3D
  equals 2D within 8.7e-15 (peak Λ̂_n 7.7). The swapped 2D MOSFET example's I_d needs the solvers: slice 3.
- QV1: the whole `nitcad_assemble_test` (179 cases), `nitcad_app_test` (77), solve "quantum confinement",
  "density gradient", interface and trap cases (21), deck confinement cases (4): all pass. Owed before the merge,
  with the owner's word: NG1 (traces byte-identical), NG2 decks, and the full Release and Debug suites.
- The 1D MOS in inversion, with only a substrate contact, has a smallest pivot ratio of 2.3e-15 without the model
  (9.8e-15 with it), below the 1e-11 floor. The QV4 Newton check on it runs with min_pivot_ratio 0.

**Slice 3a (2026-11-19, owner "choose recommended": the recommended next step, slice 3, in three local commits —
3a solver handoff, 3b transient, small signal and traces, 3c results files, plot, app, checkpoints and run-record
warnings).** As built:
- `solve_equilibrium` with the model on: the classical equilibrium first, then the stride-3 system from Λ̂ = 0
  (QP8). Every result (fields, bands, gate and trap charges) is read from the drift-diffusion state the equilibrium
  slaves (`state_from_quantum_potential`), so no equilibrium-system query has to know about Λ.
- QP9's ladder (src/solve/density_gradient.cpp). If the first solve fails with non_convergence, it runs
  k = 1/8, 1/4, 1/2, 1, each rung from the previous solution. A failed rung halves the gap to the last converged
  k, at most 4 times per rung, and then the start fails with "density gradient: the coefficient continuation failed
  at k = … (converged: …)". β·k is built as γ·k, since β is linear in γ and the barrier term is β B. The rungs
  are in `EquilibriumResult::quantum_continuation`; the run-record entry comes with 3c.
- Bias sweeps. Drift-diffusion starts from the model's equilibrium, Λ included. Quasi-static sweeps run on the
  stride-3 system, with each point's results read from the slaved drift-diffusion state at its bias. An initial
  state carries Λ (`NodeFields::quantum_potential_n_eV`, `_p_eV`, empty without the model).
- QP12, first part: the drift-diffusion E_Fn gains +Λ̂_n and E_Fp gains −Λ̂_p, and the Λ fields are in every bias
  point and the equilibrium result. They are observed only when non-empty, so every other run's trace is
  unchanged (QV1).
- Still refused until 3b: transient, small signal and traces. Adaptation is refused until its field transfer
  carries Λ.

Results (Release, targeted):
- QV4, solve part (met). The meshed MOS-C (traps and surface recombination) and the lumped one, both statistics:
  - equilibrium E_F within 1e-9 eV;
  - at −2 V and 1.5 V, drift-diffusion from the quasi-static state: the substrate current within its resolution;
  - E_F of both states within 1e-9 eV;
  - SRH and interface rates within 1e-9 of their scales (as QC5).
- The 2D MOSFET example at V_G = 1 V, equilibrium: every terminal current within its resolution and the E_F spread
  1.0e-16 eV (0.4 s, 12 iterations). Peak Λ_n is only 0.3 meV and the largest |Λ| 1.5 meV: the example's
  silicon mesh does not resolve the inversion layer, so the model has almost nothing to act on (QP15; QV13 and QV15
  report it). The case is hidden like MLDA's ([.dg]) and run by name.
- QV9, example part (met): the reflected 2D MOSFET example at V_D = 0.1 V and V_G = 0.5, 1, 1.5 V gives I_d within
  2.2e-16 of the unreflected.
- QP9: on the meshed MOS-C, Newton capped at 5 iterations fails from Λ̂ = 0 (6 suffice) and converges through the
  ladder to the uncapped state within 1e-6. Capped at 2, it fails at k = 1/128 with "converged: none".
- QV1: the whole `nitcad_assemble_test` and `nitcad_app_test`, solve confinement, density-gradient, interface,
  trap and equilibrium cases (40), deck confinement cases (5): all pass.

**Slice 3b (2026-11-20, same "choose recommended").** As built:
- Time step (QP11). Λ rows carry no storage. With incomplete ionization the dopant storage reads Λ (η_c + Λ̂_n,
  η_v + Λ̂_p, as ruled in slice 2), so the step's Jacobian gains the (n, Λ_n) and (p, Λ_p) entries, and so does
  `small_signal_matrix`. The time-step and small-signal preconditions no longer exclude the model.
- Small-signal traps. `charge_small_signal` and `InterfaceEdges::small_signal` read Λ through
  `statistics(node, x)`. They return the complex Λ partials of the trap charge and rates, and of the interface
  potential's implicit solve (Λ moves the charge only), into the semiconductor node's Λ columns, as `drift`'s.
- Starts. Bias sweeps and transients start from the equilibrium's drift-diffusion state with Λ, through one
  helper, `detail::state_from_equilibrium`. Small signal's operating points come from `sweep_bias`, which carries
  Λ.
- Fields. Transient snapshots carry Λ, and small-signal fields gain `quantum_potential_n_eV` and `_p_eV`
  (complex, eV/V; empty without the model).
- Error estimate. The transient's local error estimate measures Λ̂ absolutely in kT, as ψ, on the nodes where it
  measures ψ.
- Refusals. Transient, small signal, traces and the adjoint no longer refuse the model.
- Adaptation (ruled, recommended): it is refused at `adapt()`'s entry, `density_gradient_adaptation_refused`,
  with "quantum_confinement density_gradient: mesh adaptation does not transfer the quantum potentials". The
  field transfer does not carry Λ, and adaptation is outside Unit 62's slice plan, so this is left as a follow-up
  unit. The old "not built yet" refusal is gone.

Results (Release, targeted):
- Time-step Jacobian against central differences (house gate 1e-6), at rates 1 and 100 × t0/0.35 ps. Met on:
  meshed 1D with traps, meshed 2D with traps, meshed 2D with incomplete ionization, lumped 1D; both statistics.
- `small_signal_matrix` at s = rate equals the step's Jacobian, and at s = iω equals J + iωC without traps. Both
  within 1e-12 of the largest entry.
- QV10 (met):
  - Lumped-gate 2D MOSFET at V_G = 1.5 V, V_D = 0.1 V: Y(0)'s drain row against fourth-order DC differences
    (dV = 1e-3 V) gives g_m within 3.4e-9 and g_d within 5.4e-10.
  - MOS-C stepped to −1.5 V: the 1 µs transient equals DC within 1e-9 in ψ and Λ (peak Λ_p 0.114 eV) and 1e-7 in
    p.
- Recorded, not gated: drift-diffusion straight from the model's equilibrium to −1.5 V on the MOS-C fails
  (83 iterations, then a backward error). The DC reference there starts from the quasi-static state. QV13 reports
  convergence on the examples.
- Traces with the model have no case of their own yet; their start goes through the same `sweep_bias` path.
- QV1: the whole `nitcad_assemble_test` (180 cases) and `nitcad_solve_test` (246 passed, 4 skipped) pass, and
  every Release target builds with zero warnings.

**Slice 3c (2026-11-20, same "choose recommended").** As built:
- Files. Λ_n and Λ_p are optional keys: `NodeFields` keys 7 and 8, small-signal field keys 9 and 10 (complex). They
  are written only when the model is on, so classical files are unchanged and the format version stays (as Unit
  52's temperatures). An unpaired Λ is refused as corrupt. Steady and resume checkpoints carry Λ through the same
  fields.
- Compare. `compare_sweeps` compares Λ at the potential's tolerance (eV against V), and a Λ field missing on one
  side is a structure mismatch.
- Plot and app (QP12). New `render::Field`s: `Ln`, `Lp`, `Ec_eff` = E_c + Λ_n and `Ev_eff` = E_v − Λ_p. A result
  without the model gives an error that names the deck key. The app's field list offers them only when the result
  has Λ (as Unit 52's T_n and T_p). Ruled, recommended: the effective edges are derived in `field_values` from the
  result's bands and Λ, not stored in derived.hpp. The result already holds both, so nothing new is stored.
- Restart (QV12). A restart whose checkpoint carries Λ when the run does not solve the model, or the reverse, is
  refused in all four analyses with "the restart checkpoint '…' was solved with[out] quantum_confinement
  density_gradient and this run does[ not] solve it: start from it with initial instead". Transient and trace
  resumes are also refused by their identity. A steady restart has no identity check (Unit 61), so this is the
  check there. `initial` from a classical state stays allowed (Λ starts at 0).
- Run-record warnings. New `RunRecord::warnings` (run-record key 6, repeated text; not in the identity; printed by
  `nitcad run` as `warning: …`):
  - QP9: "the equilibrium start used the coefficient continuation (k = …)", from sweeps (drift-diffusion and
    quasi-static) and transients; small signal copies its operating points'.
  - QP15: `CoarseScan` over every result's Λ. A bulk node is a semiconductor node with no edge to an insulator
    node and not a lumped-gate node. Counted when |Λ̂| > 10, with the node count, the results, the largest value
    and its node. Traces do not scan yet.
- QV12 vocabulary pins: the deck's `[models]` list and the app's completion list grew in slice 1 (24 completion
  items) and are unchanged here.

Results (Release, targeted):
- QV12, lumped-gate MOS-C deck (0.5 nm surface spacing), V_G 0 to −1.2 V in −0.1 V steps (accumulation; peak Λ_p
  0.117 eV):
  - the sweep file round-trips Λ bit for bit;
  - `Ln`, `Ec_eff` and `Ev_eff` equal the result's Λ and bands;
  - restarted from the −0.6 V checkpoint, the remaining 6 points match the uninterrupted sweep: ψ within 2.6e-16,
    Λ within 1e-17 eV, currents within their resolution (a MOS-C at DC);
  - the model change is refused both ways.
- Not met as worded: QV12 asks for a sweep restart that is "bit-identical … (as Unit 61 RG)". Unit 61 made only
  transient and trace resumes bit for bit. Its steady restart is to 1e-9 (DK9, RG6), because a steady checkpoint
  holds physical-unit fields, not the solver state. Measured here at 1e-15 to 1e-16, not to the bit. The criterion
  is frozen, so it is put to the owner by name: rule QV12's sweep restart to Unit 61's steady standard (1e-9), or
  ask for a bit-for-bit steady resume (a new resume record).
  **Ruled (2026-11-20, owner "choose recommended"): (a).** QV12's sweep restart is held to Unit 61's steady
  standard, 1e-9 relative (DK9, RG6), which is met at 1e-16. The transient resume with the model on, bit for bit,
  is the "as Unit 61 RG" case. QV12 is met.
- The same "as Unit 61 RG" check that is bit for bit, a transient resume with the model on: the deck diode pulse,
  stopped after 3 and after n/2 steps and resumed, equals the uninterrupted run bit for bit (times, steps, orders,
  error ratios, currents and Λ in the last snapshot). The plain deck is refused on the same checkpoint. The
  diode's Λ is below 1e-17 eV, so this case checks the resume path, not Λ's size.
- Recorded: on this MOS-C, drift-diffusion in inversion stops at V_G = 0.4 V classically too (36 iterations at
  0.3 V, then 100). This is the 1D MOS-C's minority supply, not the model; density gradient stops one step
  earlier.
- Warnings: the texts, bulk selection (a lumped-gate node at 20 kT not counted) and count pinned by a unit case. A
  clean sweep records none. Results-file round-trip of Λ and the warnings: `nitcad_io_test` passes.
- QV1:
  - pass: `nitcad_assemble_test`, `nitcad_app_test`, `nitcad_render_test`, `nitcad_analysis_test`,
    `nitcad_io_test`;
  - `nitcad_solve_test`: 247 passed, 4 skipped;
  - `nitcad_deck_test`: 115 passed, plus its 5 existing failed-as-expected Schottky cases;
  - every Release target builds with zero warnings.

**Slice 4, part 1 (2026-11-20, owner "choose recommended and start slice 4").** QV6–QV8 are in
tests/verify/benchmark/density_gradient_test.cpp, hidden ([.dg_validation], Release, run by name; 1 s in all).
The setup:
- A 1D meshed MOS-C: n+ poly, built-in SiO2 (the reference's barrier 3.0 eV and mass 0.5), p-type, silicon to
  0.3 µm, 20 oxide cells, the interface straddled at h = 0.1 nm, spacing h to 3 nm, then ×1.08 up to 2 nm.
- The quasi-static sweep, so E_F = 0, as in the reference.
- Boltzmann statistics, and bandgap narrowing off (the reference has none). V_FB comes from the classical gate
  charge's zero.
- The reference has the same material constants as NiTCAD's silicon, valleys (0.916, 0.19) ×2 and
  (0.19, √(0.916·0.19)) ×4, and its slice-0 grid.
- N_inv = ∫ n over silicon boxes, and the centroid = ∫ x n / N_inv.
- Ruled (recommended): the reference's guard N_A < 0.01 N_v becomes N_A < N_v, because QV6 asks for 1e18. With
  Boltzmann statistics on both sides the bulk only needs E_F in the gap. The reference's own tests pass
  unchanged (6 cases).

Results:
- Classical cross-check (gated): NiTCAD against the reference's classical electrons gives V_th within 0.2 mV and
  N_inv within 0.35% in all 4 devices. With bandgap narrowing on, 1e18 differed by 27–34 mV, which is how the
  BGN setting was found.
- QV6, reported (no verdict until QP0). γ = 1 and MLDA's masses (placeholders). Errors are against
  Schrödinger–Poisson, at V_th(classical SP) + 0.5, 1, 1.5 V:

  | t_ox, N_A | N_inv DG (MLDA) | centroid DG (MLDA) [nm] | ΔV_th DG, MLDA, SP [mV] |
  |---|---|---|---|
  | 1.5 nm, 1e17 | +11.9, +9.2, +8.0% (+9.1, +8.6, +8.4%) | −0.57, −0.46, −0.40 (≈ −0.48) | 20.0, 29.1, 46.1 |
  | 1.5 nm, 1e18 | +17.0, +11.2, +9.1% (+14.3, +10.9, +9.9%) | −0.45, −0.39, −0.36 (≈ −0.48) | 45.2, 58.9, 97.8 |
  | 5 nm, 1e17 | +5.8, +4.5, +3.8% (+3.9, +3.3, +3.0%) | −0.64, −0.57, −0.52 (≈ −0.47) | 21.0, 30.6, 48.5 |
  | 5 nm, 1e18 | +11.5, +7.4, +5.7% (+8.6, +5.8, +4.7%) | −0.48, −0.45, −0.42 (≈ −0.47) | 54.7, 71.8, 118.8 |

  Against the placeholder bounds (10%, 0.3 nm, 25%), DG would miss:
  - the centroid at every point (its centroid lies 0.36–0.64 nm nearer the surface);
  - ΔV_th everywhere, at 54–57% of SP's;
  - N_inv at 4 of 12 points.

  DG's confinement here is about half the reference's in ΔV_th, and MLDA does better on all three except the thin
  1e18 centroid. This is what QP0's coefficient (Ancona's r factor) has to settle: no tuning, as QV6 says. The
  mass mismatch (one isotropic 0.98 against the per-valley masses) and Λ_p's lack of a reference stand as stated.
- QV7 (met): the 5 nm, 1e17 device at V_G − V_FB = 2 V gives N_inv 3.084834e12, 3.082945e12 and 3.082551e12
  at 0.2, 0.1 and 0.05 nm. The last two differ by 1.3e-4 (≤ 1%), less than the first two (6.1e-4).
- QV8 (met):
  - x_p → 0 (m_ox 0.5 → 128, x_p 0.159 → 0.010 nm): the surface density falls monotonically, 3.30e18 → 3.45e17,
    and the peak moves inward, 0.7 → 0.9 nm.
  - γ → 0: the difference from classical, gated over γ = 1e-3, 1e-4, 1e-5, has order 0.992 and 0.999 in N_inv
    and 0.968 and 0.997 in max |n − n_cl|.
  - The pair from 1e-2 is reported, not gated: order 0.949 and 0.810. There the surface layer still moves by
    20%, so it is outside the limit's range. The bound is unchanged; only the γ sequence was extended to 1e-5.

**Slice 4, part 2: the examples (2026-11-20).** QV13 and QV15 are in tests/deck/quantum_confinement_test.cpp
"the examples on their meshes, and refined" ([.dg_examples], Release, run by name, about 14 min). Each example
deck's own sweep, classical, MLDA and DG. V_th is the largest-g_m extrapolation minus V_D/2; the swing is the
smallest ΔV/Δlog I_d.
- QV15 (met on the examples' own meshes: every point converges):
  - 1D diode: largest |Λ̂| 0.023 at 0.005 µm, I(0.7 V) −7.0e-5, 73 against 72 iterations.
  - 2D MOSFET: 0.093 at (1.5, 0.5) µm (the drain junction), I_on −4e-4.
  - 3D MOSFET: 2.93 at (1.4, 0.5, 0.4) µm (the drain junction), I_on −2e-4.
  - No QP15 warning: nothing reaches 10 kT away from the interfaces.
  - On these meshes (10 nm at the surface in 2D) the model does almost nothing. The largest Λ sits at the
    junctions, not in the channel.
- The 2D MOSFET refined to 0.1 nm at the surface (10 nm to 0.94 µm, 1 nm to 0.948 µm, 0.1 nm to 0.95 µm), reported:
  - Largest |Λ̂| 3.54 at the interface (1.45, 0.9499) µm.
  - With the deck's defaults (no step cutting), DG stops after 14 of 16 points (V_G = 1.4 V: Newton not converged
    in 100 iterations). Classical and MLDA converge.
  - With `max_cuts = 8` (a deck option, the same for every model) DG converges at every point:

    | refined 2D MOSFET | V_th [V] | swing [mV/dec] | I_on [A/cm] |
    |---|---|---|---|
    | classical | 0.5031 | 95.8 | 7.7314e-2 |
    | MLDA | 0.5100 | 96.0 | 7.6521e-2 |
    | DG | 0.5075 | 95.9 | 7.6811e-2 (−0.65%) |

  - The refined mesh shows the coarse runs' near-zero effect was under-resolution (QP15), not the model being
    absent. DG's ΔV_th (+4.4 mV) is again below MLDA's (+6.9 mV), consistent with QV6.
- QV13 (reported):

  | example | classical | MLDA | DG |
  |---|---|---|---|
  | 2D MOSFET: V_th, swing, I_on | 0.5707 V, 99.9, 6.5274e-2 | the same to 7 digits | 0.5712 V, 99.9, 6.5251e-2 |
  | 2D MOSFET: time, iterations | 0.75 s, 83 | 0.73 s, 83 | 2.39 s, 97 |
  | 3D MOSFET: V_th, swing, I_on | 0.4849 V, 107.9, 6.9463e-6 | the same to 7 digits | 0.4849 V, 107.9, 6.9446e-6 |
  | 3D MOSFET: time, iterations | 4.44 s, 39 | 4.42 s, 39 | 38.6 s, 80 |
  | refined 2D MOSFET, `max_cuts = 8` | 2.36 s, 100 | 2.48 s, 102 | 744 s, 242 accepted (QP9 continuation not used) |

  - MLDA equals classical on the coarse examples: their nearest silicon nodes are far from the oxide (Unit 47
    QC12).
  - Cost per accepted iteration against drift-diffusion: 2.7× in 2D, 4.2× in 3D.
  - The refined DG time includes every failed attempt of the cut steps (100 iterations each), which the
    iteration count does not show. On this mesh, DG's Newton needs smaller bias steps than classical: a cost
    finding, not a fix.

**QP0 recorded (2026-11-20, owner supplied the papers, "(a)").** The four sources were read from the owner's
copies, outside the repo. In our form the quantum potential is Λ = −(γħ²/6m) ∇²√n/√n. In Ancona's form it is
2b ∇²√n/√n (in volts). So γ/m = 12bq/ħ².

- **(i) Coefficient.**
  - Ancona and Tiersten, PRB 35 (1987) 7959, eq. (4.9a) and p. 7964:
    - the static DG equation is a "Schrödinger equation" with coefficient 2b;
    - matching the one-electron effective-mass ground state gives b = ħ²/(4qm*), with m* the normal mass;
    - the "note added in proof" says a more rigorous microscopic derivation gives b = ħ²/(12m*q).
  - Ancona and Iafrate, PRB 39 (1989) 9536, eq. (3.12): b = ħ²/(12m*q). This is the lowest-order Wigner
    correction for a Boltzmann, high-temperature gas of infinite extent; eq. (3.13) gives the condition
    ε = ħ²β/(8m*d²) ≪ 1. In our form it is **γ = 1** with m = m*.
  - Ancona, IEEE TED 47 (2000) 1449, eq. (2.7): b = ħ²/(4qm*r), where r = 1 is the pure state (one subband)
    and r = 3 the high-temperature limit, which is the 1989 value.
  - Ancona 2000, pp. 1452–1453 and Fig. 7: for silicon inversion layers m* and r vary oppositely, and he
    takes m*r = 1.
  - Wettstein, Schenk and Fichtner, IEEE TED 48 (2001) 279, eq. (11), has the same Λ form as ours with γ as
    a fit factor: 3.6, fitted on p. 281. It is not a theoretical value, so it is not the default (TVM2).
  - The placeholder (γ = 1 with m = 0.98) paired the high-temperature coefficient with the pure-state normal
    mass. That made it about 3× too weak: γ/m = 1.02 where Ancona's m*r = 1 gives 3. This is QV6's earlier
    under-confinement.
- **(ii) Masses.**
  - Electrons: Ancona 2000, p. 1451, gives the single-valley bulk average (m_t² m_l)^{1/3} = 0.328. Footnote 1:
    m*r = 0.328 × 3 = 0.984, his recommended m*r = 1.
  - **Default: γ = 1 (1989 eq. (3.12)) with m_n = 0.328.** That is γ/m = 3.05, within 2% of m*r = 1.
  - Holes: none of the four sources gives a value. m_p stays MLDA's 0.49, marked as unsourced. Λ_p has no
    reference anyway (QP14).
  - Wettstein et al. use the density-of-states mass (p. 281). Not adopted, because their γ is fitted with it.
- **(iii) Oxide boundary condition.**
  - Ancona and Tiersten 1987, eq. (4.4b), and Ancona 2000, p. 1451: ρ = 0 at the surface, an insulating oxide
    and an infinite barrier, with tunnelling omitted.
  - Wettstein et al.: the barrier enters through eq. (8), with ∇(Φ+Λ)·n = 0 at the outer boundaries.
  - **QP5 confirmed (a), not revised.** Our reference (QP14) has a finite barrier, so the penetration form
    matches the physics QV6 compares against. QV8 shows the hard wall is its x_p → 0 limit (QR3 (c) is that
    limit).
- **(iv) Published agreement, which sets QV6's bounds.**
  - Ancona 2000, p. 1454, against SHRED Schrödinger–Poisson (16 subbands, Fermi–Dirac; 31 Å oxide; N_A 1e15,
    3.8e16 and 1e18):
    - V_th within 10 mV, worst at the heaviest doping;
    - N_inv "quite accurate" and the centroid close, shown in figures only, with no number in the text.
  - Wettstein et al. publish agreement only for their fitted γ, so it does not apply to the theoretical
    coefficient.
  - **Bounds frozen:**
    - |(ΔV_th DG) − (ΔV_th SP)| ≤ 10 mV, cited;
    - |ΔN_inv/N_inv| ≤ 10% and |Δcentroid| ≤ 0.3 nm, the slice-0 numbers kept as engineering bounds. The
      sources give no number for these two, and this is stated.
    - The 25% relative shift placeholder is replaced by the cited 10 mV.

**QV6 verdict (QP0's defaults; tests/verify/benchmark/density_gradient_test.cpp, now gated).**

| t_ox, N_A | shift DG / SP [mV] | miss [mV] | N_inv DG vs SP | centroid DG − SP |
|---|---|---|---|---|
| 1.5 nm, 1e17 | 35.9 / 46.1 | 10.2 | +2.0 to +3.6% | −0.09 to −0.16 nm |
| 1.5 nm, 1e18 | 80.2 / 97.8 | 17.6 | +2.3 to +4.8% | −0.06 to −0.09 nm |
| 5 nm, 1e17 | 37.6 / 48.5 | 10.9 | +1.2 to +2.1% | −0.14 to −0.20 nm |
| 5 nm, 1e18 | 97.0 / 118.8 | 21.8 | +1.7 to +3.7% | −0.08 to −0.10 nm |

- **N_inv and the centroid pass at all 12 points. The threshold shift fails at all four devices**: DG gives
  78–82% of SP's shift, missing by 10.2–21.8 mV against the 10 mV bound.
- The placeholder γ gave 54–57% of SP's shift, with N_inv up to +17% and the centroid up to −0.64 nm.
- MLDA (unchanged): 60–63% of SP's shift.
- Recorded as a fail, not tuned (QV6).
- Known differences from Ancona's comparison, stated rather than corrected:
  - our reference is Boltzmann, with a finite barrier and two valley sets (m_l = 0.916);
  - his is Fermi–Dirac, an infinite barrier and m_l = 0.98;
  - V_th here is N_inv = 1e12 cm⁻².
- Ancona's Fig. 7 (m*r = 0.9 "slightly better") suggests a stronger coefficient would narrow the gap. Choosing
  one is calibration, so it is **the owner's ruling** (QV6: "the owner may then rule a reported γ"). Not done.
- QV7 under the new defaults: 3.2e-4 and 7.3e-5 (pass).
- QV8 under the new defaults: orders 0.899/0.980/0.998 (N_inv) and 0.669/0.919/0.990 (n). The heavier
  coefficient makes the 1e-2 pair more pre-asymptotic; it was already report-only. Pass.
- Code changes:
  - physics/quantum_confinement.hpp: `dg_default_electron_mass` = 0.328, with the citations;
  - the app's deck vocabulary and sections: the `dg_mass_n` default is 0.328;
  - the physics and deck default pins were updated.
- The slice-4 part 2 example numbers (QV13, QV15) and the FinFET runs were measured with the placeholder
  mass. They are rerun under the new default before the unit closes.

**Slice 4 examples rerun under QP0's default (2026-11-20, m_n = 0.328; tests/deck/quantum_confinement_test.cpp
[.dg_examples], all 26 assertions pass).**

| example | classical | MLDA | DG |
|---|---|---|---|
| 2D MOSFET: V_th, swing, I_on | 0.5707 V, 99.9, 6.5274e-2 | the same | 0.5720 V, 99.9, 6.5205e-2 (I_on −0.11%) |
| 2D MOSFET: time, iterations | 0.69 s, 83 | 0.68 s, 83 | 2.17 s, 96 |
| 3D MOSFET: V_th, swing, I_on | 0.4849 V, 107.9, 6.9463e-6 | the same | 0.4850 V, 107.9, 6.9413e-6 (I_on −0.07%) |
| 3D MOSFET: time, iterations | 4.17 s, 39 | 4.14 s, 39 | 29.3 s, 66 |
| refined 2D, `max_cuts = 8`: V_th, swing, I_on | 0.5031 V, 95.8, 7.7314e-2 | 0.5100 V, 96.0, 7.6521e-2 | 0.5114 V, 96.1, 7.6386e-2 (I_on −1.20%) |
| refined 2D, `max_cuts = 8`: time, iterations | 2.14 s, 100 | 2.21 s, 102 | 910 s, 283 |

- QV15 is met on the examples' own meshes again. The largest |Λ| is 0.12 kT (2D) and 4.3 kT (3D), both
  below the 10 kT warning.
- On the refined 2D mesh, DG still stops at point 5 with the default `max_cuts`. With 8 it completes, the
  largest |Λ| 4.2 kT.
- The stronger, cited coefficient moves the refined DG V_th from +4.4 mV to +8.3 mV above classical
  (+1.4 mV above MLDA).
- QV13: DG costs 2.7× per iteration in 2D and 4.2× in 3D, unchanged. The refined case's 910 s again includes
  the failed cut attempts.
- **The FinFET is dropped (owner, 2026-11-20: "dont run on finfet, run on mosfet").**
  - The classical and MLDA FinFET runs completed: 446 s and 397 s, on the placeholder build.
  - The DG run was first refused, correctly (QP5): the deck's `insulator` spacer has no band data.
  - With band data added to the run deck (Si₃N₄ and HfO₂, Robertson 2000), it ran 19 min without finishing
    point 1 and was stopped. That is a cost finding, not a verdict.
  - The FinFET is not gated (QP14). The example deck's spacer would need band data for DG (QP10 says so).

**Ruling (2026-11-20, owner "choose recommended for γ for QV6"): QR4 (a) stands.**
- The default stays the cited theoretical coefficient, γ = 1 with m_n = 0.328. No calibrated default is
  shipped (TVM2).
- QV6 stays recorded as a **fail** on the threshold shift: DG gives 78–82% of SP's shift, missing by
  10.2–21.8 mV against the 10 mV bound. It passes on N_inv and the centroid.
- The gated test keeps its 10 mV check, so `[.dg_validation]` reports the miss whenever it is run. It is
  hidden and not in ctest.
- A deck may set `dg_gamma_n` / `dg_mass_n`, and the run record always reports them. That is the place for a
  user's calibration, not the default.
- The bound is not loosened. Loosening it would be a new owner ruling, by name.

### Unit 58 — solver scale-up: AMG, coupled preconditioning, deterministic threads, parallel assembly (proposed 2026-11-15; on `linalg/solver-scale-up`)

**Owner request (2026-11-15):** "Plan for solver scale-up (AMG, parallel assembly)". This is the 11.3 row "Scalable
linear solvers", next step "Algebraic multigrid and block (field-split) preconditioners; parallel assembly". The
branch starts from `main` `beb494f`. The number 58 was used only informally for the MoS₂ gap list of 2026-11-12,
which was never recorded here.

**What is measured today (the starting point):**
- **GMRES-ILUT (Unit 44)**, one Newton step's linear solve, one thread: 3D diode 24k / 81k / 192k unknowns took
  0.77 / 3.8 / 14.1 s, with **40 / 60 / 90 iterations**. The growth with size is what ILUT cannot remove. SparseLU
  took 1021 s and 10.9 GB at 192k.
- **PARDISO (Unit 18)** at 81k 3D unknowns: 4.7 s on one thread, 1.07 s on eight. It is the only multithreaded part.
  The FinFET deck (45,375 nodes) takes 10.6 min on 6 threads (Unit 50).
- **Assembly is below 1% of a Newton iteration in 2D and 3D** (Unit 52, energy-balance measurement); the
  factorization dominates. Once a faster preconditioner shrinks the solve, that share rises, but it is unmeasured.
- GMRES-ILUT requires `threads = 1`. Its products, dot products and triangular solves are serial.
- GMRES fails the complex small-signal systems (LQ3a: resolution up to 93% of max |Y|), which stay on SparseLU.

**Standing conditions:** opt-in. The default backend and every existing result stay untouched. The 6.8 rules hold:
the same input, thread count and build give the same result, and a threaded kernel parallelizes over the output
index and never accumulates in thread order. The result is also independent of the core it runs on (Unit 41).
`linalg` gains no dependency on `assemble`, and no Eigen type enters its public headers (Q5's check). Criteria are
fixed before code. Release and Debug suites pass at every commit with `/W4 /WX`, plus hand mutants. Nothing is
pushed or merged without authorization.

**Left out of this unit:** GPU backends (cuDSS was dropped 2026-11-15); MPI and distributed memory; MUMPS; complex
(small-signal) AMG, since LQ3a stays; making any new path the default (SQ6); a natural-test Newton norm; dropping the
insulator nodes' unknowns (FD5).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| SD1 | 0 | **Measure first.** A Release `[.benchmark]` profile of one Newton iteration, split into assembly, analysis, factorization or preconditioner setup, solve, update and the rest. The cases: the 2D MOSFET, the 3D diode at 30³ and 40³, the reduced swept FinFET, the tensor FinFET deck, and a 3D energy-balance device. Each is run with SparseLU (where it fits), PARDISO at 1 and 6 threads, and GMRES-ILUT. A 3D diode size series then runs up to about 1e6 unknowns with GMRES-ILUT and PARDISO, recording iterations, time and memory. Measured pinned to P cores; recorded here before slice 1. | No solver change. This sets the numbers SG4–SG7 are judged against. |
| SD2 | 1 | **Deterministic threads.** A persistent pool in `base` (`std::jthread`) with `parallel_for` over **fixed-size chunks, independent of the thread count**. Reductions are per-chunk partial sums, combined in chunk order, so the bits are identical for any thread count. Used for the CSR product, GMRES's dot products, norms and updates, the equilibration scaling, and the backward-error check. `threads > 1` is accepted for `gmres_ilut`. The ILUT triangular solves stay serial. | 6.8's rule, by construction. The gate-tunnelling workers (GR11) may move onto the pool later; not in this unit. |
| SD3 | 2 | **Algebraic multigrid, in-tree (SQ1).** Smoothed aggregation (Vaněk, Mandel, Brezina 1996): a strength-of-connection threshold θ, serial greedy aggregation in node order (deterministic), a damped-Jacobi smoothed prolongator, and Galerkin coarse operators. Damped Jacobi takes l1 diagonal scaling, the smoother is Chebyshev (both parallel over rows, so deterministic), and the coarse level is solved by SparseLU below a size limit. It is a V-cycle used as GMRES's right preconditioner. Aggregates are built once per pattern (`analyze`), and the operators are rebuilt per Jacobian (`factorize`). It runs on the equilibrated matrix, like ILUT. | The same in-tree choice as GMRES (Unit 44): no new dependency, and determinism under our control. |
| SD4 | 3 | **Coupled preconditioning (SQ2): CPR-style two-stage.** Stage 1 is AMG on the potential block, A_ψψ, with a reduction of the carrier columns into it. Stage 2 is ILUT on the full system for the residual left by stage 1. The block labels reach `linalg` as plain data (`std::span<const std::int32_t>`, one field label per unknown); the solve layer gives them from `DofMap` (Unit 40). Every equation set works (3, 4 or 5 unknowns per node), and insulator nodes are labelled ψ. | The potential equation is elliptic, which AMG handles well; the continuity rows are convection-dominated, which ILUT handles. AMG on the whole system is measured as SD3's own backend for comparison. |
| SD5 | 4 | **Parallel assembly (SQ4), deterministic, conditional.** The terms compute their per-edge and per-node contributions in parallel into per-item buffers. A row-parallel gather then adds them into F and J **in the serial order of every entry** (FE3), so the result is bit-identical to serial assembly for every thread count, and FD6's objection does not apply. Built only if SD1 or the post-SD4 profile shows assembly at **≥ 10% of a Newton iteration** on the 3D cases, and only after Unit 57 (terms then read `SystemData` with no shared mutable state). Otherwise the measurement is recorded and the slice dropped. | The buffers cost about 48 doubles per edge for drift-diffusion; reported. |
| SD6 | — | **Configuration and deck (SQ7).** `SolverConfig::preconditioner = ilut \| amg \| cpr` under a `gmres` backend; `gmres_ilut` stays accepted and means `gmres` + `ilut` with an unchanged run identity. The new options are AMG's θ, smoother degree, coarse size and maximum levels. The deck adds `[solver] preconditioner = …`, `amg_*` keys and `threads` for `gmres`. The run record digests the new settings only when used. Small signal stays on SparseLU (LQ3a). | |
| SD7 | — | **Commits.** One branch. Commit 0 is SD1, then one commit per slice. The owner may stop after any slice; slice 1 is useful alone (threaded GMRES-ILUT). ARCHITECTURE.md 6.10 and 11.3 are updated to "as built" after acceptance. | |

**Criteria (fixed before code; a target missed is recorded as a limitation, never loosened):**

| # | Criterion |
|---|---|
| SG1 | **Default untouched:** with the default backend, and with `gmres_ilut` at one thread, every suite and the FG1/FG2 traces are unchanged bit for bit. |
| SG2 | **Determinism:** for SD2–SD5's kernels and solves, results are bit-identical at 1, 2, 4, 6 and 8 threads, and across P and E cores (a pinned probe as `linalg_eigen_blocking`); repeats are bit-identical. |
| SG3 | **Accuracy:** each new preconditioner meets `max_backward_error` on the original system (the layer's own check, as Unit 44). Device cases (diode, MOS C-V, transient, 2D MOSFET, a 3D deck, electrothermal, energy balance) agree with SparseLU or PARDISO within Unit 44's `[gmres]` bounds. Fallback and `fallback = none` behave as for ILUT. |
| SG4 | **Scalability:** on the 3D diode series, AMG on the equilibrium Poisson systems and CPR on the bias Jacobians grow GMRES iterations by **at most 1.5×** from 24k to the largest size in SD1 (ILUT measured 2.25× from 24k to 192k). |
| SG5 | **Speed at one thread:** on the 3D diode at 40³ and SD1's largest case, the per-Newton-step linear time with CPR is **≤ 0.5×** GMRES-ILUT's. |
| SG6 | **Speed with threads:** at 6 threads, the linear time of SG5's cases is **≥ 2.5×** faster than at one thread, and the tensor FinFET deck's sweep (17 points) is faster with `gmres` + `cpr` at 6 threads than with PARDISO at 6 threads. Its currents agree with PARDISO's within the deck's terminal-current resolution. |
| SG7 | **Memory:** the new preconditioners use no more peak memory than ILUT on SG5's cases. |
| SG8 | **Parallel assembly (if built):** F and J bit-identical to serial for 1–8 threads, and `evaluate` ≥ 3× faster at 6 threads on the 3D cases. |
| SG9 | **Suites and mutants:** Release and Debug suites pass at every slice, with no warnings. At least 8 mutants per slice are caught. At least: a reduction combined in thread order; a chunk size derived from the thread count; aggregation iterating in a hashed order; a missing Galerkin transpose; the smoother's damping off by a factor; CPR's stage 2 skipped; ψ labels shifted by one; the gather adding in edge-buffer order rather than entry order. |
| SG10 | **Layering:** `linalg` has no dependency on `assemble` (labels are plain integers); there is no Eigen type in public headers, and no new vcpkg dependency (SQ1 a). |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| SQ1 | AMG implementation | (a) In-tree smoothed aggregation (SD3). (b) AMGCL (header-only, MIT, in vcpkg; a new dependency whose own dependencies and threaded determinism need checking). (c) hypre BoomerAMG (C, MPI-oriented, large). | **(a)**, as Unit 44 did for GMRES: no dependency, and 6.8's determinism is ours to guarantee. |
| SQ2 | Coupled preconditioner | (a) CPR two-stage: AMG on ψ, ILUT on the full system. (b) Block Gauss–Seidel over ψ, n, p with AMG on each block. (c) Point-block AMG on the whole system (3×3 to 5×5 node blocks). | **(a)**: the standard for coupled elliptic and convection-dominated systems; it reuses ILUT; (c) is measured as SD3's whole-system AMG. |
| SQ3 | Threading mechanism | (a) The in-tree fixed-chunk pool (SD2). (b) OpenMP (MSVC's OpenMP 2.0; reduction order not guaranteed). | **(a)**: the only way 6.8's "identical across thread counts" holds by construction. |
| SQ4 | Parallel assembly | (a) Conditional: built only if assembly is ≥ 10% of an iteration after SD4, and after Unit 57. (b) Build it unconditionally. (c) Drop it. | **(a)**: measured at below 1% today; (b) risks a slice with no gain. |
| SQ5 | Order against Unit 57 | (a) Slices 0–3 now (they need only `DofMap`, already public); slice 4 after Unit 57. (b) Everything after Unit 57. | **(a)**: no dependency, and slices 0–3 touch only `linalg`, `solve` and `deck`, not the assembler Unit 57 rewrites. |
| SQ6 | Default | (a) Opt-in; a default change is a separate owner ruling after the measurements. (b) `gmres` + `cpr` becomes the default for 3D. | **(a)**: as for every solver unit so far. |
| SQ7 | Deck syntax | (a) `backend = gmres` with `preconditioner = ilut \| amg \| cpr`; `gmres_ilut` kept as an alias with an unchanged identity. (b) One backend name per combination (`gmres_amg`, `gmres_cpr`). | **(a)**: options compose without new names, and old decks are unchanged. |

**Risks:** AMG on badly scaled drift-diffusion rows is the known hard case; CPR exists because plain AMG on the
coupled system often fails, so SG4–SG6 may be missed and recorded as limitations. The FinFET target (SG6) competes
with a mature multithreaded direct solver at only about 136k unknowns. The gain grows with size, so SD1's size series
decides whether this pays off on today's decks or only on larger ones.

**Rulings (2026-11-20, owner "Unit 58 choose recommended"):** SQ1 (a) in-tree smoothed aggregation, no new
dependency; SQ2 (a) CPR two-stage, AMG on ψ and ILUT on the full system; SQ3 (a) the in-tree fixed-chunk pool;
SQ4 (a) parallel assembly only if assembly is ≥ 10% of an iteration after SD4; SQ5 (a) slices 0–3 first, slice 4
after Unit 57 (Unit 57 merged on `main` 2306228, so that condition is met); SQ6 (a) opt-in, any default change a
separate owner ruling after the measurements; SQ7 (a) `backend = gmres` with `preconditioner = ilut | amg | cpr`,
`gmres_ilut` kept as an alias with an unchanged run identity. The proposal was rebased onto `main` 3c36224 (it was
written on beb494f). Not started: slice 0 waits for the owner's word.

**Slice 0 (SD1), 2026-11-20 (owner "start slice 0"). As built:** `tests/deck/solver_profile_test.cpp` in
`nitcad_deck_test` (it reads decks), hidden: `[.scale_profile]` (the phase split) and `[.scale_series]` (the size
series), Release, run pinned to the six P cores (`start /b /wait /affinity 3C3`, CPUs 0, 1, 6–9 of the Core Ultra 5
225F; 32 GB). No product code changed. Each case starts from its bias stamped on its equilibrium state, as
`sweep_bias`'s first point, and is measured twice per configuration:
- the **split**: a plain Newton loop in the test (full damped steps, no line search) timing assembly (`evaluate`),
  the analysis, the first and later factorizations (ILUT builds for GMRES-ILUT), the solves, the update and the rest;
- **`newton_solve`** itself, with the case's own Newton options, on a forwarding system that times `evaluate` and
  the update; the remainder is the linear solver plus the loop's O(n) work.
GMRES-ILUT runs strict (`fallback = none`) so a failure shows instead of hiding in SparseLU. SparseLU only where it
fits (2D, and the 30³ diode capped at 3 iterations).

Results, per Newton iteration (split; "fact." = factorization or ILUT build):

| Case (unknowns) | Backend | Its. | Assembly | Fact. | Solve | Assembly share (split / `newton_solve`) |
|---|---|---|---|---|---|---|
| 2D MOSFET example (3,198) | GMRES-ILUT | 5 | 0.1 ms | 6.0 ms | 3.5 ms (27 products) | 1.4% / 1.2% |
| | PARDISO 1 / 6 thr. | 6 | 0.1 ms | 2.9 / 0.9 ms | 0.2 / 0.1 ms | 2.7% / 2.5%; 4.8% / 4.6% |
| | SparseLU | 6 | 0.1 ms | 7.6 ms | 0.4 ms | 1.5% |
| 3D diode 30³ (81,000) | GMRES-ILUT | 10 | 5.5 ms | 2.38 s | 0.72 s (56) | 0.2% / 0.2% |
| | PARDISO 1 / 6 thr. | 10 | 6.0 ms | 4.82 / 1.10 s | 0.04 / 0.02 s | 0.1% / 0.1%; 0.5% / 0.5% |
| | SparseLU (3 its.) | 3 | 5.7 ms | 77.4 s | 0.14 s | 0.0% (held 3.1 GB) |
| 3D diode 40³ (192,000) | GMRES-ILUT | 10 | 14.6 ms | 6.53 s | 2.56 s (79) | 0.2% / 0.2% |
| | PARDISO 1 / 6 thr. | 10 / 11 | 15 ms | 27.7 / 7.42 s | 0.11 / 0.07 s | 0.1% / 0.1%; 0.2% / 0.2% |
| 3D energy balance, EC19's diode on 41×21×21 at −0.4 V (90,405, stride 5) | GMRES-ILUT | — | 32 ms | 5.68 s | fails | GMRES does not reach the backward error in 1000 products at Newton iteration 2 |
| | PARDISO 1 / 6 thr. | 3 | 14–17 ms | 3.91 / 0.87 s | 0.04 / 0.02 s | 0.3% / 0.3%; 1.7% / 1.7% |
| Reduced swept FinFET (24,900) | GMRES-ILUT | 6 | 2.0 ms | 0.12 s | 0.03 s (25) | 1.2% / 1.2% |
| | PARDISO 1 / 6 thr. | 7 | 2.0 ms | 0.38 / 0.09 s | 7 / 5 ms | 0.5% / 0.5%; 2.0% / 2.1% |
| Tensor FinFET deck, first point (136,125) | GMRES-ILUT | 21 | 13 ms | 1.64 s | 1.22 s (79) | 0.5% / 0.8% |
| | PARDISO 1 / 6 thr. | 27 / 23 | 13 ms | 10.6 / 2.72 s | 0.07 / 0.03 s | 0.1% / 0.4%; 0.5% / 4.2% |

- **Assembly is at most 2% of a Newton iteration in every 3D case** (4.8% on the tiny 2D example with PARDISO at 6
  threads). The `newton_solve` figure of 4.2% (FinFET deck, PARDISO 6 threads) comes from its line search: 533
  evaluations over 56 iterations from this start, against 13 iterations and 25 evaluations with GMRES-ILUT and 48
  evaluations with PARDISO at one thread. The paths differ only in rounding, so this start is near the edge of the
  deck's line search (the deck itself steps the gates in 0.05 V and converges). Recorded, not investigated here.
- **GMRES-ILUT's time is the ILUT build (57–76%) and the solves (20–42%)**, both serial; the analysis is 0.2%.
  PARDISO's is its factorization (93–99.8%).
- SD5's condition (assembly ≥ 10% of an iteration) is not met before SD4; it is measured again after SD4.

Size series (one Newton step's linear solve on LG4's 3D diode at 0.5 V, a fresh solver: analysis, factorization or
ILUT, one solve; peak memory as the solver's held private bytes):

| Size (unknowns) | GMRES-ILUT: time (ILUT), products, held | PARDISO 1 thread | PARDISO 6 threads: time, held |
|---|---|---|---|
| 20³ (24,000) | 0.64 s (0.53), 40, 47 MB | 0.46 s, 169 MB | 0.18 s, 254 MB |
| 30³ (81,000) | 3.17 s (2.42), 60, 156 MB | 5.21 s, 806 MB | 1.49 s, 842 MB |
| 40³ (192,000) | 9.51 s (6.61), 90, 376 MB | 29.3 s, 2.7 GB | 8.40 s, 2.7 GB |
| 50³ (375,000) | 21.2 s (14.1), 108, 735 MB | — | 42.6 s, 7.1 GB |
| 60³ (648,000) | 39.1 s (25.6), 116, 1.3 GB | — | 138 s, 15.3 GB |
| 70³ (1,029,000) | 68.4 s (42.0), 142, 2.0 GB | — | not measured: the run was stopped by the host for low memory during this factorization (fill 1.78e9 at 60³) |

- GMRES iterations grow **3.55×** from 24k to 1.03M unknowns (40 → 142): the growth SG4 (≤ 1.5×) is judged against.
- PARDISO at 6 threads is faster up to 40³ and slower from 50³ (2× at 50³, 3.5× at 60³); its memory grows about as
  n^1.6 and passes 32 GB before 1e6 unknowns on this machine. GMRES-ILUT's grows linearly (2 GB at 1e6).
- SG5 and SG7's references: the 40³ GMRES-ILUT step is 9.51 s, 376 MB; the largest (70³) 68.4 s, 2.0 GB.

**Slice 1 (SD2), 2026-11-20 (owner "After slice 0 is done, choose recommended and start slice 1").**
- **Ruling under "choose recommended": reductions stay serial.** SD2 as written combines per-chunk partial sums of
  the dot products and norms in chunk order. That order is not the one-thread loop's, so `gmres_ilut` at one thread
  would change in the last bits and SG1 would fail. Computing them serially at one thread only would make one thread
  differ from two, and SG2 would fail. The frozen criteria win over the design text: the pool runs only per-index
  kernels (each output entry computed in its own fixed order), the dot products and norms stay serial, and every
  thread count gives the one-thread bits, which are today's. Chunked reductions remain possible later only with a
  new owner ruling on SG1's one-thread bits.
- **As built.** `base::ThreadPool` (`include/NiTCAD/base/thread_pool.hpp`): `threads − 1` persistent `std::jthread`
  workers plus the caller; `parallel_for(n, chunk, body)` over fixed chunks `[k chunk, min(n, (k+1) chunk))`, claimed
  through an atomic counter; workers spin for 2 ms after a job, then sleep on a condition variable (the first version
  woke them through the condition variable every time, which made GMRES slower at 6 threads than at 1: 3.92 s against
  2.97 s at 40³). `base::for_chunks(pool, n, chunk, body)` is one call over `[0, n)` without a pool.
  `gmres_ilut` accepts 1 to 1024 threads (was exactly 1); with more than one the solver owns a pool, and GMRES runs on
  it the matrix product (chunks of 2048 rows), the basis normalizations, the Gram–Schmidt updates, w = Σ tᵢ vᵢ (each
  entry summed over i in order) and the iterate update (chunks of 8192 entries). Serial: the ILUT build and its
  triangular solves, the dot products and norms, the backward-error check and the equilibration. The deck's
  `[solver] threads` already reaches `SolverConfig::threads`; the `preconditioner` key comes with SD6.
- **SG1 (one thread unchanged).** The new test's one-thread solve (160² convection–diffusion, restart 8, 63
  products) gives the same FNV-1a digest of x, iterations and backward error as the code before the change, built
  from `main` `e57918a` (whose `linalg` is this branch's base): `5c5ca5bc96ed4618`, 63, `0x1.e769a3bc79725p-29`.
  The full suites and FG1/FG2 traces are not run (owner's word needed).
- **SG2 (determinism).** `gmres_ilut: every thread count gives the one-thread bits`: x, iterations and backward error
  identical at 2, 3, 4, 6, 8 threads (25,600 unknowns, several chunks and restarts); the complex case at 2, 5, 8
  threads. `[.scale_threads]` (40³ and 60³ diode steps): every entry identical at 1 and 6 threads. The pool's own
  tests: every index exactly once and the fixed chunks for 1, 2, 3, 6, 8 threads and chunk sizes 1, 64, 8192, reused
  across calls; chunks run on more than one thread. The cross-core probe (P against E cores) is not run: the pooled
  kernels have no cache-dependent blocking (Unit 41's cause was Eigen's), so the bits cannot depend on the core.
- **Speed (`[.scale_threads]`, P cores): no useful gain.** GMRES solve 40³: 2.98 s at 1 thread, 3.06 s at 6; 60³:
  13.80 s, 12.98 s. The ILUT build (6.7 s and 26 s) is unchanged. What the pool can run is a small part of a GMRES-ILUT
  solve: the ILUT triangular solves (about 5× the matrix's nonzeros) and the serial dot products dominate, and the
  vector updates are memory-bound. The pool is the base of slices 2–3, whose smoothers and products are per-index.
- **SG9 (mutants): 9 of 9 caught** (script `C:\Users\disha\bisect-tools\mutants58s1.sh`, M4 and M8 by hand): the
  last partial chunk dropped; chunks starting one late; the claim counter not reset between jobs; workers never
  started (the caller runs every chunk); the product skipping each chunk's first row; each update chunk dropping its
  last entry; w's sum starting at i = 1; `gmres_ilut` refusing more than one thread; the pool running everything on
  the caller. SG9's reduction mutants (thread-order combination, chunk size from the thread count) do not apply: the
  pool has no reduction, and a per-index kernel gives the same bits for any chunk size. A solver that ignores its
  pool gives the same bits too; only the timing shows it.
- **Tests run (Release, targeted):** `nitcad_base_test` (15 cases), `nitcad_linalg_test` (83), `nitcad_solve_test
  "[gmres]"` (4), `[.scale_threads]`.
- **Full suites (owner "run full suite and debug tests", 2026-11-20), at `b4342aa`:** Release 251/251 pass (808 s);
  Debug 230/230 pass (2251 s). The Debug build first failed on `main`'s own code, not this branch's:
  `tests/solve/quantum_confinement_test.cpp` (Unit 62, merged as PR #102 before its CI finished) had a `FAIL` before
  a `return` in a lambda, which is C4702 (unreachable code) under `/WX` in Debug. Fixed here (the lambda returns an
  error whose message the test's check rejects) and the Debug suite ran with the fix; the Release run predates the
  fix, which only touches that test. On the rebase onto main 737cfc3 (2026-11-20) this file takes main's own fix
  (PR #103, 9fdffca: `FAIL_CHECK` before the return), so the branch no longer changes it. Not run: NG1-style traces.

**Slice 2 (SD3), 2026-11-20 (owner "start slice 2 of Unit 58"; standing "choose recommended").**
- **As built.** `src/linalg/amg.{hpp,cpp}` (private, Eigen-free): smoothed-aggregation AMG as GMRES's right
  preconditioner, one V-cycle per application. `SolverConfig::preconditioner = ilut | amg` (default `ilut`) with
  `amg_strength_threshold` (θ, default 0.08, [0, 1)), `amg_smoother_degree` (2, 1–10), `amg_coarse_size` (1000,
  1–100000), `amg_max_levels` (10, 2–50), validated for every backend. `amg` needs the `gmres_ilut` backend and a
  real matrix (complex is invalid_input: small signal stays on SparseLU, LQ3a). The run record digests
  `linear.preconditioner` and the `amg_*` keys only when the preconditioner is not ILUT (SG1: identities unchanged).
  The deck keys (SD6) are not in this slice. A hierarchy that cannot be built (a coarsest level SparseLU rejects, a
  non-finite entry) is inaccurate_solve, so the fallback and `fallback = none` behave as for ILUT. Setup is serial;
  the V-cycle's products and vector updates run on the slice-1 pool, per index in fixed chunks.
- **Design rulings under "choose recommended"** (SD3's text against measurements; no criterion touched):
  1. *Similarity scaling.* SD3 said "runs on the equilibrated matrix, like ILUT". R A C is not symmetric when A is, and
     Galerkin with R = Pᵀ then failed: a heterogeneous Poisson problem (conductances 2⁻¹⁰–2¹⁰) did not converge in
     1000 iterations on R A C. AMG is built on the similar matrix U (R A C) U⁻¹ = G A G (U = √(C/R), G = √(RC)) and
     applied as U⁻¹ V(U v); the near-null vector is 1/G (the constant of the original unknowns).
  2. *Weighted lumping.* Weak entries are lumped into the filtered diagonal weighted by the near-null vector
     (a_ii += a_ij b_j / b_i), so A_F b = A b. With plain lumping the equilibrated heterogeneous case took 146
     iterations against 20 unscaled; weighted, 22.
  3. *Diagonal, not l1.* SD3's l1 scaling is not invariant under the symmetric scaling (120 against 46 iterations).
     The prolongator smoother and the Chebyshev smoother use the true diagonal; λ_max(D⁻¹A) is estimated by 15 power
     iterations from a fixed pseudo-random vector (serial, so bit-reproducible), times 1.1; ω = 4/(3λ).
  4. *Chebyshev interval [λ/10, λ]* (not λ/30): 2D Poisson 129² 10 iterations against 14, heterogeneous 17 against 22.
  5. *A stalled coarsening is smoothed, not factorized.* On the equilibrium Poisson systems (diagonally dominant) the
     first version made every weakly coupled node a singleton, the coarsening stopped near the full size, and SparseLU
     got that level (40³: 7.0 M stored against 0.8 M). A level that stalls (no aggregate, or shrinking by less than
     10%) above `coarse_size` is now only smoothed. Leaving unconnected nodes out of aggregates instead (Vaněk's
     isolated nodes) measured worse on coarse levels (48³ Poisson 16 against 12 iterations) and was not kept.
- **Tests (`tests/linalg/amg_test.cpp`, 13 cases):** transpose, product and the pooled product against dense
  arithmetic; the coarsening of a chain (aggregates, Dirichlet rows isolated, repeatable); strength symmetry and θ;
  Rᵀ = P and A_c = PᵀAP by dense arithmetic on every level; the configuration checks; GMRES-AMG strict on 2D
  Poisson 33²–257² (9, 9, 10, 10 iterations) and 3D 16³–48³ (9, 10, 12), each checked by an independent backward
  error; the equilibrated heterogeneous case (17; homogeneous 10); convection 0.5 (14); degrees 1–5 and level limits;
  a stalled diagonally dominant matrix (one smoothed level, nothing stored); one analysis per pattern, one hierarchy
  per matrix, repeat bits; the coarsening kept for a later matrix of the same pattern; a singular coarsest level
  (strict: inaccurate_solve; fallback: SparseLU's singular_system); 1, 2, 3, 4, 6, 8 threads bit-identical (200²,
  16 iterations, several cycles). Linalg suite (Release): 96 cases pass.
- **Limitation (recorded, not hidden): whole-system AMG fails on convection-dominated and coupled systems.** Upwind
  convection at cell Péclet 2 and the coupled badly scaled fixture do not converge in 1000 iterations (the fallback
  solves them), nor does the 3D diode's drift-diffusion Jacobian at 0.5 V at 20³, 30³ or 40³ (below). This is the case
  CPR (SD4, slice 3) exists for: AMG on the ψ block only.
- **Measured (`[.scale_amg]`, Release, pinned to the P cores; one fresh solve: build and solve).** The equilibrium
  Poisson Jacobian of the 3D diode (one unknown per node), at the charge-neutral guess and at the converged potential
  (unit right-hand side), strict:

| Size (unknowns) | GMRES-ILUT, first step: total, its, held | GMRES-AMG, first step | GMRES-ILUT, converged | GMRES-AMG, converged: 1 / 6 threads |
|---|---|---|---|---|
| 20³ (8,000) | 0.02 s, 4, 9 MB | 0.01 s, 7, 5 MB | 0.02 s, 3 | 0.01 s, 8 / 0.01 s |
| 30³ (27,000) | 0.08 s, 5, 30 MB | 0.05 s, 12, 24 MB | 0.09 s, 4 | 0.05 s, 9 / 0.04 s |
| 40³ (64,000) | 0.26 s, 5, 70 MB | 0.11 s, 11, 35 MB | 0.29 s, 6 | 0.13 s, 9 / 0.09 s |
| 50³ (125,000) | 0.63 s, 6, 136 MB | 0.25 s, 11, 73 MB | 0.70 s, 7 | 0.25 s, 9 / 0.19 s |
| 60³ (216,000) | 1.57 s, 7, 236 MB | 0.59 s, 12, 121 MB | 1.71 s, 9 | 0.53 s, 9 / 0.40 s |
| 70³ (343,000) | 2.45 s, 7, 382 MB | 0.87 s, 12, 182 MB | 3.10 s, 10 | 0.87 s, 9 / 3.90 s (see below) |

  - **SG4 (AMG part), equilibrium Poisson:** iterations grow 7 → 12 (**1.71×, misses ≤ 1.5×**) at the first Newton
    step and 8 → 9 (**1.13×, passes**) at the converged potential, from 8,000 to 343,000 unknowns. The first step's
    growth is the 20³ point (7); from 30³ on it is 12, 11, 11, 12, 12. ILUT on these systems grows 4 → 7 and 3 → 10:
    the equilibrium Poisson rows are strongly diagonally dominant, so they are easy for both. SG4's CPR part is slice 3.
  - **Speed and memory on these systems:** AMG is 2.6–3.6× faster than ILUT from 40³ and holds half the memory.
  - **Threads:** 6 threads gave 1.3–1.5× on the AMG solve up to 60³ (bits identical, checked in the unit test). At 70³
    the 6-thread solve took 3.35 s against 0.41 s at one thread: a single run, not repeated; recorded as observed
    and to be re-measured with slice 3's thread runs before any SG6 judgement.
  - **Drift-diffusion at 0.5 V, whole-system AMG** (SD3's comparison): fails strict at 20³, 30³, 40³ (1000 products;
    2.8, 12.0, 35.9 s); GMRES-ILUT solves them in 40, 60, 90 products (0.88, 4.16, 10.87 s), as in slice 0.
- **SG9 (mutants): 17 of 17 caught** (scripts `C:\Users\disha\bisect-tools\mutants58s2.sh` and `mutants58s2b.sh`): SG9's
  aggregation in a hashed order, a missing Galerkin transpose (R from the tentative P), the smoother's and the
  prolongator's damping off by 2; and the equilibration not passed to AMG, unweighted lumping, the coarse correction
  dropped, no post-smoothing, λ underestimated by 2, the coarsening rebuilt per matrix, Dirichlet rows with explicit
  zeros aggregated (survived the first run; caught after a test with explicit-zero contact rows was added), one-way
  strength, a coarse LU failure ignored, the product overwriting instead of adding, AMG never selected, complex AMG
  accepted, a stalled level factorized. SG9's CPR and gather mutants belong to slices 3 and 4.
- **Not run:** the full Release and Debug suites (owner's word needed; Debug builds also need it), FG1/FG2 traces,
  the cross-core probe (as slice 1: no cache-dependent blocking in the AMG kernels).

**Slice 3 (SD4), 2026-11-20 (owner "start slice 3"; standing "choose recommended").**
- **As built.** `src/linalg/cpr.{hpp,cpp}` (private, Eigen-free): `Preconditioner::cpr` under `gmres_ilut`, real
  only (complex is invalid_input), its two stages taking the `amg_*` and `ilut_*` settings. One application is
  z1 = P AMG(R v), z = z1 + ILUT(v − A z1), on the equilibrated matrix as ILUT: R injects the potential rows, P
  prolongs a potential correction to each node's other unknowns. The labels reach `linalg` as plain integers:
  `BasicLinearSolver::set_field_labels(std::span<const std::int32_t>)` (0 potential, positive a field coupled to it
  at the same node, negative an unknown in no node). A node is a potential and the positive labels right after it.
  Labels that differ from the last ones forget the analysis; the same labels keep it. Without one label per unknown,
  a cpr factorization is invalid_input. `assemble::field_labels(const DofMap&)` derives them. `newton_solve` sets
  them when the system has `dofs()` and the preconditioner is cpr, so equilibrium, bias, transient and trace Newton
  runs get them. Also wired: the trace's bordered solver (the arc unknown labelled −1), `adapt_adjoint`'s two solvers,
  and the weighted-current Newton solver. The weighted-current Laplace solver takes `amg` instead (one unknown per
  node, no fields). The run record already digests `linear.preconditioner` (2 for cpr) and the `amg_*` keys only
  when the preconditioner is not ILUT. Deck keys (SD6) are not in this slice.
- **Fix of slice 2's code:** small signal copied the steady config and switched only the backend to SparseLU, so with
  `amg` (and now `cpr`) it failed at create ("a preconditioner other than ilut needs the gmres backend"). It now
  resets the preconditioner to ILUT too. The bug was latent: no deck key reached it yet.
- **Design ruling under "choose recommended": the reduction.** SD4 said "AMG on A_ψψ with a reduction of the carrier
  columns into it". A_ψψ alone, in the (ψ, n, p) unknowns, is λ² times the Laplacian: it lacks the screening term,
  which dominates where the device is doped. Left (quasi-IMPES) decoupling with the carrier diagonals does not
  restore it either: at equilibrium the ψ-derivatives of a Scharfetter–Gummel continuity row sum to zero, so it adds
  a zero-row-sum operator. Built instead: A_p = R A P, with P's entries the node-local response
  e = −B_cc⁻¹ B_cs (B the node's own block). For Scharfetter–Gummel at equilibrium this is exactly e_n = n and
  e_p = −p (an identity of the Bernoulli function), so A_p is the nonlinear Poisson Jacobian, λ² Laplacian +
  (n + p), with A_ψψ's pattern. It is generic for 3–5 unknowns per node. A node whose B_cc has an exactly zero pivot
  (or gives a non-finite e) is left uncoupled (e = 0). Insulator and contact nodes give e = 0 by themselves (their
  carrier rows hold the densities). R A P with the same P as the prolongation means an exact stage 1 removes the
  potential rows' residual before ILUT. The AMG gets the potential unknowns' equilibration scales, since the
  equilibrated A_p is R_ψ A_p C_ψ.
- **Tests (`tests/linalg/cpr_test.cpp`, 9 cases)** on an independent 2D Scharfetter–Gummel (ψ, n, p) Jacobian
  (pn junction, doping 10⁴ nᵢ, λ² = 0.05, contacts as identity rows, optional quasi-Fermi bias):
  - the responses at equilibrium are n and −p (to 1e-12), and every entry of A_p is λ² L + (n + p), contacts identity;
  - the configuration and label checks: SparseLU + cpr, complex, no labels, short labels, a coupled label before any
    potential, no potential; labels ignored by ILUT;
  - GMRES-CPR strict on 24², 48², 96² at equilibrium and biased, each against an independent backward error, with
    fewer iterations than ILUT where ILUT converges;
  - a bordering unknown labelled −1; the same labels keeping the analysis and new ones redoing it; the fallback and
    strict failure; a singular node block giving e = 0;
  - 1, 2, 3, 4, 6, 8 threads bit-identical (12,288 unknowns, 53 iterations).
  - given equilibration scales decades apart (rows and columns), the AMG's finest level is the symmetric G A_p G
    where A_p is symmetric (permittivities 2^-10 to 2^10). Iteration counts cannot show this: stage 1's accuracy
    barely moves them (below).

  Solve layer (`[gmres]`): a new bias case (labels from the system, one analysis over two Newton runs, J(0.5 V) and
  the 0–0.7 V sweep against SparseLU: relative 1e-10, worst 0.000 of max(1e-10 |J|, resolution)); the MOS C-V,
  transient switch-on and small-signal junction cases rerun with cpr: gate charge within 7.9e-17, currents within
  3.5e-21 of the largest, admittance 0.117 of its bound. `fields: the cpr field labels`. Linalg suite: 105 cases pass.
- **The fixture shows where the limit is.** GMRES-CPR / GMRES-ILUT iterations at equilibrium: 24, 60, 330 against
  32, 690, failed at 24², 48², 96²; biased: 14, 36, 120 against 17, 111, failed. CPR beats ILUT, but it still grows
  about 2–5× per halving of h. A probe with an exact stage 1 (the coarsest level the whole of A_p) gave the same
  counts (327 against 330 at 96²); a tighter ILUT (drop 1e-6, fill 40) gave 41. **The growth is stage 2's.** The
  continuity blocks are diffusion operators, elliptic too, and a one-level ILUT does not scale on them. CPR's premise,
  one elliptic field and the rest local or hyperbolic, holds for reservoir pressure, not for drift-diffusion
  carriers.
- **Measured (`[.scale_cpr]`, Release, pinned to the P cores):** one Newton step's linear solve on the 3D diode at
  0.5 V as `[.scale_series]` (fresh solver, strict):

| Size (unknowns) | GMRES-ILUT: total (build), products, held | GMRES-CPR 1 thread | GMRES-CPR 6 threads |
|---|---|---|---|
| 20³ (24,000) | 0.64 s (0.53), 40, 47 MB | 0.67 s (0.53), 42, 55 MB | 0.72 s, 42 |
| 30³ (81,000) | 3.21 s (2.44), 60, 156 MB | 3.40 s (2.50), 60, 189 MB | 3.56 s, 60 |
| 40³ (192,000) | 9.70 s (6.68), 90, 376 MB | 9.78 s (6.83), 71, 452 MB | 9.33 s, 71 |
| 50³ (375,000) | 21.7 s (14.3), 108, 735 MB | 21.9 s (14.5), 96, 883 MB | 21.2 s, 96 |
| 60³ (648,000) | 39.7 s (25.9), 116, 1.27 GB | 41.7 s (26.2), 112, 1.53 GB | 39.5 s, 112 |
| 70³ (1,029,000) | 68.5 s (42.1), 142, 2.02 GB | 72.4 s (42.6), 138, 2.42 GB | 72.7 s, 138 |

  - **SG4 (CPR part): missed.** 42 → 138 products, **3.3×** from 24k to 1.03M (ILUT 3.55×; target ≤ 1.5×).
  - **SG5: missed.** 40³: 9.78 s against 9.70 s (**1.01×**, target ≤ 0.5×); 70³: 72.4 s against 68.5 s (1.06×). The
    time is ILUT's build (stage 2 builds the same ILUT as plain GMRES-ILUT) plus solves, which CPR does not shorten.
  - **SG6 (linear part): missed.** 6 threads: 1.05× at 40³, 1.00× at 70³ (target ≥ 2.5×); the serial ILUT build and
    triangular solves dominate. The FinFET deck sweep needs the deck keys (SD6) and was not run; with these
    numbers it would not beat PARDISO at 6 threads.
  - **SG7: missed.** CPR holds ILUT's memory plus A_p and its hierarchy: +17–20% (452 against 376 MB at 40³).
  - **SG2 (threads): met.** Bits identical in the unit test at 2–8 threads; the 6-thread runs here gave the same
    products and backward errors.
  - **SG3: met** on the device cases above.
  - The 70³ AMG six-thread anomaly of slice 2 (3.90 s against 0.87 s, one run) **did not repeat**: three repeats gave
    0.15–0.16 s solve at 6 threads against 0.33–0.37 s at one (2.2×), bits as at one thread. Recorded as a one-off.
- **What would address it (not built; each a new owner ruling, since SQ2 ruled (a)):** multigrid for the carrier
  blocks too, that is SQ2 (b) (block Gauss–Seidel over ψ, n, p with AMG on each block, the continuity blocks
  nonsymmetric) or (c) (point-block AMG on 3×3 node blocks). Plain AMG on the whole system failed in slice 2. A
  tighter ILUT alone moves the counts (41 against 330 on the fixture) but not the growth or the serial build.
- **SG9 (mutants): 15 of 15 caught** (script `C:\Users\disha\bisect-tools\mutants58s3.sh`): SG9's CPR stage 2
  skipped and ψ labels shifted by one (in `field_labels`), and the prolongation into the next node; the responses'
  sign flipped, the responses dropped (A_ψψ alone), the prolongation of the potential only, stage 2's residual not
  updated, the scales not given to the AMG (survived the first run; caught after the symmetry test was added), new
  labels keeping the analysis, a singular node block kept (first written so it did not compile), complex cpr
  accepted, cpr never selected, a coupled label before its potential accepted, `newton_solve` giving no labels, small
  signal keeping cpr. SG9's gather mutant belongs to slice 4.
- **Not run:** the full Release and Debug suites (owner's word needed), FG1/FG2 traces, the FinFET deck sweep
  (SG6, needs SD6's deck keys).

**Wrap-up, 2026-11-20 (owner "choose recommended" after slice 3, then "go with item 1").** The recommendation
chosen: build SD6's deck keys, record slice 4 as not built under SQ4 (a), do not start multigrid on the carrier
blocks (SQ2 b or c would be a new proposal with its own criteria fixed before code), and update ARCHITECTURE.md.
- **SD6 (deck), as built.** `[solver] backend = gmres` is `gmres_ilut` (SQ7 a: the same backend and run identity;
  `gmres_ilut` stays accepted). `preconditioner = ilut | amg | cpr` (default `ilut`); `amg_strength_threshold`
  ([0, 1)), `amg_smoother_degree` (1–10), `amg_coarse_size` (1–100000), `amg_max_levels` (2–50); `threads` already
  reached `gmres`. A preconditioner other than ilut with another backend is a deck error. The app's deck vocabulary
  (`backends`, `preconditioners`, five solver keys) and the section menu carry the new keys. Tests: `deck build:
  backend = gmres, the preconditioner and the amg keys` (gmres with an unchanged identity, ilut explicit unchanged,
  amg and cpr selected, each amg key in the identity with cpr, the range and backend errors), the Unit 44 key test's
  unknown-backend message, the app's deck sections and vocabulary tests (solver keys 25 → 30). Targeted runs pass:
  those two deck-build cases, `app deck sections*` (4), `app vocabulary*` (8).
- **Slice 4 (SD5, parallel assembly): not built.** SQ4 (a) builds it only if assembly is ≥ 10% of a Newton iteration
  after SD4. Slice 0 measured it at ≤ 2% on every 3D case with GMRES-ILUT, and CPR does not shorten the linear part
  (1.01–1.06× GMRES-ILUT's time per step, slice 3), so the share stays at ≤ 2%. No separate profile was run.
- **SG6's FinFET deck sweep:** not run. The deck keys now allow it, but slice 3's measurements already show that
  GMRES-CPR has no thread gain, so it cannot beat PARDISO at 6 threads.
- **ARCHITECTURE.md:** 6.10 gains "As built (Unit 58)" (the pool, the preconditioners and their deck keys, AMG, CPR
  and its measured limit, parallel assembly not built). The 11.3 row now says what Unit 58 built; its next step is
  multigrid on the carrier blocks.
- **Unit 58's open results:** SG1–SG3 met; SG4–SG7 missed (recorded, not loosened); SG8 not applicable (no parallel
  assembly); SG9 met per slice (mutants 9, 17, 15); SG10 met (labels are plain integers, no Eigen in public headers,
  no new dependency). The full Release and Debug suites have not run on slices 2–3 and this wrap-up (owner's word
  needed).

### Unit 63 — multigrid on the carrier blocks (Unit 58 follow-up, SQ2 b) (proposed 2026-11-20; on `linalg/carrier-multigrid`)

**Owner request (2026-11-20):** "start Unit 58 follow-up, choose recommended". This is the next step Unit 58's wrap-up
recorded for the 11.3 row "Scalable linear solvers": multigrid for the carrier blocks, SQ2 (b) or (c), as a new
proposal with its own criteria fixed before code. The branch starts from `main` `2065afc`.

**What Unit 58 measured (the starting point):**
- GMRES-CPR (AMG on the reduced potential system, then ILUT of the full matrix) grows 42 → 138 products from 24k to
  1.03M unknowns on the 3D diode at 0.5 V (3.3×; GMRES-ILUT 3.55×). It is no faster than GMRES-ILUT and has no
  thread gain, since ILUT's serial build (57–76% of the time) and triangular solves dominate.
- On the 2D fixture (`cpr_test.cpp`: 24², 48², 96²) GMRES-CPR takes 24/60/330 iterations at equilibrium and 14/36/120
  biased. An exact stage 1 gives the same counts: **the growth is stage 2's.** The continuity blocks are diffusion
  operators, elliptic too, and a one-level ILUT does not scale on them.
- Plain AMG on the whole coupled system fails (slice 2).

**Why the carrier blocks can take the existing AMG.** A Scharfetter–Gummel continuity block is not symmetric, but it
is a diagonal similarity away from a symmetric one: its entries satisfy a_ij / a_ji = c_i / c_j (c = e^{ψ/V_T} for
electrons, e^{−ψ/V_T} for holes, up to the rows' and columns' scales), so with D = diag(c)^{-1/2} the matrix
S = D A D⁻¹ is symmetric. Its off-diagonal entries are −g Δ / (2 sinh(Δ/2)), g the edge's coefficient and Δ the edge's potential step in V_T:
bounded, exponentially weak across high-field edges (which become diagonally dominant rows, which smoothing handles).
The near-null vector of S is 1/d (the equilibrium density, transformed). D can be computed from the matrix alone,
with no potential, so `linalg` still sees only values and labels.

**Standing conditions:** as Unit 58: opt-in; the default backend and every existing preconditioner (`ilut`, `amg`,
`cpr`) are untouched; 6.8's determinism (the same bits for every thread count); `linalg` gains no dependency on
`assemble` and no Eigen type in a public header; criteria fixed before code; Release and Debug suites with
`/W4 /WX` before a push (the suites wait for the owner's word); hand mutants; nothing pushed or merged without
authorization.

**Left out:** point-block AMG (SQ2 c); a nonsymmetric AMG (AIR, Petrov–Galerkin); complex (small-signal) AMG;
parallel assembly (SD5 stays not built); making any new path the default; the setup's parallelization.

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| MD1 | 0 | **Measure the carrier blocks alone first.** A test-only probe: the electron and hole blocks of the Unit 58 fixture (24², 48², 96²; bias 0 and 20) and of the 3D diode series' 0.5 V Jacobians (20³–70³, equilibrated as GMRES sees them), each solved by GMRES with (i) the slice-2 AMG on the block as it is and (ii) AMG on the symmetrized block (MD2). Iterations, setup and solve time, levels. Recorded here before slice 1. | No product code beyond MD2. Go/no-go: if (ii) does not meet MG5's 1.5× growth on the blocks alone, slice 1 is not started and the owner is asked. |
| MD2 | 0 | **Algebraic symmetrization** (`linalg/detail`, with the AMG). For a square block A: ln d_j − ln d_i = ½ ln(a_ij / a_ji) along a spanning forest of its off-diagonal pairs (breadth-first in node order from each component's lowest unknown; a pair is used only if both entries are non-zero and of the same sign). Each component's mean ln d is shifted to 0, and ln d is clamped to ±230 (d within 1e±100) so reverse-biased blocks do not overflow. S = D A D⁻¹ and the hierarchy is built on its symmetric part ½(S + Sᵀ) (equal to S when exact). Given to the AMG as matrix S with scales R = C = d, which makes its similarity the identity and its near-null vector 1/d. One application: z = D⁻¹ V(D r). | Exact for Scharfetter–Gummel blocks (any mobility, recombination on the diagonal); approximate where impact ionization or Fermi–Dirac terms break the ratio, which only weakens the preconditioner. A block with no off-diagonal pairs gets d = 1. |
| MD3 | 1 | **Preconditioner `cpr_amg`.** Stage 1 as `cpr` (AMG on the reduced potential system, prolonged by the local responses). Stage 2, on the residual stage 1 leaves: one forward block Gauss–Seidel sweep over the field labels in order (0, 1, 2, …): z_f = AMG_f⁻¹(r_f − Σ_{g<f} A_fg z_g), each AMG_f one V-cycle on the symmetrized diagonal block of label f (MD2). Unknowns outside every node (label −1) get a Jacobi step. No ILUT (MQ3). | The ψ block's AMG in stage 2 removes stage 1's inexactness on the potential rows. Setup serial; the block products and V-cycles row-parallel on the pool, as Unit 58. |
| MD4 | 1 | **Configuration and deck.** `Preconditioner::cpr_amg`; `[solver] preconditioner = cpr_amg` under `backend = gmres`. The `amg_*` keys apply to every hierarchy. The run record names it as the other preconditioners. Small signal stays on SparseLU (LQ3a). | |
| MD5 | 2 | **Measure.** MG5–MG8 on the fixture, the 3D diode series (`[.scale_cpr_amg]`, as `[.scale_cpr]`), and the tensor FinFET deck's 17-point sweep against PARDISO at 6 threads. | |
| MD6 | — | **Commits.** One branch; the proposal, then one commit per slice, then ARCHITECTURE 6.10 and 11.3 "as built" and the records. The owner may stop after any slice. | |

**Criteria (fixed before code; a target missed is recorded as a limitation, never loosened):**

| # | Criterion |
|---|---|
| MG1 | **Nothing else moves:** the default backend and `ilut`, `amg` and `cpr` give the same bits as on `main` (their existing tests unchanged, and the `[gmres]`, `[amg]` and `[cpr]` cases' printed counts equal). |
| MG2 | **Determinism:** `cpr_amg` bit-identical at 1, 2, 4, 6 and 8 threads, and on repeats. |
| MG3 | **Accuracy:** `max_backward_error` met on the original system; the device cases Unit 58 ran with `cpr` (the bias case against SparseLU, the MOS C-V, transient switch-on, small-signal junction) agree within the same bounds with `cpr_amg`. Fallback and `fallback = none` behave as for `cpr`. |
| MG4 | **Symmetrization:** on the fixture's carrier blocks S is symmetric to 1e-13 relative per entry, and d_j / d_i matches e^{∓Δψ/2} (an independent derivation from the fixture's potential) to 1e-12; a block with a sign-mismatched pair, a zero partner, two components and a reverse bias past the clamp each give finite d and the documented fallback. |
| MG5 | **Scalability:** GMRES iterations with `cpr_amg` grow **at most 1.5×** from 24² to 96² on the fixture (bias 0 and 20) and from 24k to 1.03M unknowns on the 3D diode series (Unit 58's SG4 target; `cpr` 3.3×). |
| MG6 | **Speed at one thread:** on the 3D diode at 40³ and 70³ the linear time per Newton step is **≤ 0.5×** GMRES-ILUT's (SG5's target). |
| MG7 | **Speed with threads:** at 6 threads MG6's cases are **≥ 2.5×** faster than at one thread, and the tensor FinFET deck's sweep is faster with `gmres` + `cpr_amg` at 6 threads than with PARDISO at 6 threads, currents within the deck's terminal-current resolution (SG6's targets). |
| MG8 | **Memory:** peak memory no more than GMRES-ILUT's on MG6's cases (SG7). |
| MG9 | **Mutants:** at least 8 per slice caught, at least: the half in ½ ln(a_ij/a_ji) dropped; the transform applied as D⁻¹ A D; the near-null vector d instead of 1/d; the clamp removed (the reverse-bias test overflows); a sign-mismatched pair used; Gauss–Seidel's lower couplings skipped (Jacobi); stage 1 skipped; a label's block given another label's hierarchy. |
| MG10 | **Layering:** as SG10: plain integer labels, no Eigen in public headers, no new dependency. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| MQ1 | Coupled scheme | (a) CPR's stage 1 kept, stage 2 a forward block Gauss–Seidel with AMG per field (MD3). (b) Block Gauss–Seidel alone, no CPR stage. (c) Point-block AMG on the whole system (SQ2 c). | **(a)**: stage 1 already carries the ψ–carrier coupling (the screening term, slice 3); (c) needs a new aggregation and block smoothers, and whole-system AMG failed. |
| MQ2 | AMG on a carrier block | (a) Symmetrize algebraically (MD2), then the slice-2 AMG. (b) A nonsymmetric AMG (AIR, Petrov–Galerkin smoothed aggregation). (c) The slice-2 AMG on the block as it is. | **(a)**: exact for Scharfetter–Gummel and reuses Unit 58's hierarchy; (c) is slice 0's baseline; (b) is a new method. |
| MQ3 | ILUT in stage 2 | (a) None: the per-field AMGs only. (b) ILUT of the full matrix after the sweep. | **(a)**: ILUT's serial build is what made `cpr` no faster; (b) is measured only if (a) misses MG5 on the fixture, and recorded. |
| MQ4 | Deck syntax | (a) `preconditioner = cpr_amg`. (b) A key `cpr_stage2 = ilut \| amg`. | **(a)**: one value, as SQ7 (a); `cpr` keeps its identity. |
| MQ5 | Fields other than n, p (carrier temperatures, density-gradient potentials) | (a) The same per-field symmetrized AMG. (b) n and p only; the others a Jacobi step. | **(a)**: labels are generic; a symmetric block gets d = 1. Reported per equation set. |
| MQ6 | Default | (a) Opt-in; a default change is a separate owner ruling. (b) Default for 3D. | **(a)**, as every solver unit. |

**Risks:** the high-field rows make the symmetrized blocks nearly diagonal, so the coarsening may stall there (slice
2's "stalled level smoothed" path); the block Gauss–Seidel ignores the carriers' coupling back into ψ beyond stage 1,
so the outer GMRES may still grow with size; the serial setup (aggregation, Galerkin products) may cap MG7, as the
setup was 1.71× AMG's first-step time in slice 2.

**Rulings (2026-11-20, owner "start Unit 58 follow-up, choose recommended"):** MQ1 (a) CPR stage 1 with a block
Gauss–Seidel stage 2 of per-field AMG; MQ2 (a) algebraic symmetrization; MQ3 (a) no ILUT; MQ4 (a)
`preconditioner = cpr_amg`; MQ5 (a) every field; MQ6 (a) opt-in. Slice 0 started.

**Slice 0 (MD1, MD2), 2026-11-20. As built:** `amg::symmetrize`, `amg::symmetric_part` and `SymmetrizedAmg` in
`src/linalg/amg.{hpp,cpp}` (not reachable from any configuration yet); `tests/linalg/block_probe.hpp` (a field's
diagonal block; GMRES(30) to a componentwise backward error of 1e-8, the slice-2 AMG at its default settings, one
thread). Tests in `cpr_test.cpp`: MG4 on the fixture's carrier blocks (8², 24²: S symmetric to 1e-13, d_l / d_k =
e^{∓Δψ/2} from the fixture's potential to 1e-12, 2m + 1 components, 2m unpaired entries), the unusable pairs, the
components and the clamp (a chain spanning ln d 600: 6 unknowns clamped, d² and 1/d² finite), GMRES with the
symmetrized AMG on the 24² blocks; hidden `[.carrier_blocks]` (fixture) and, in `solver_profile_test.cpp`,
`[.scale_carrier_blocks]` (3D diode) and `[.carrier_diag]` (the near-null choice). `nitcad_linalg_test`: 108 cases
pass (MG1: the `amg` and `cpr` cases unchanged).
- **A finding fixed during the slice:** the first `symmetric_part` averaged an entry with any stored partner. The 3D
  device Jacobian stores explicit zeros in the contact rows, so every coupling into a contact was halved: the ψ block
  (exactly symmetric, d = 1) then took 20–28 iterations against the plain AMG's 9–12, and the hole blocks failed at
  every size. Averaging only usable pairs (MD2's own definition) fixed both. Near-null choice on 20³ / 40³ after the
  fix (`[.carrier_diag]`): 1/d (MD2) 10/13 electrons, 11/17 holes; constant 10/12, 17/23; d 16/27, 60/148. MD2's
  1/d kept.
- **Fixture blocks (`[.carrier_blocks]`), iterations (24² is a single level, coarse_size 1000):**

| Block | 48² | 96² | 192² |
|---|---|---|---|
| ψ, AMG as is / symmetrized | 9 / 9 | 9 / 9 | 9 / 9 |
| electrons, AMG as is / symmetrized | 81 / **11** | 294 / **13** | failed / **13** |
| holes, AMG as is / symmetrized | 146 / **15** | failed / **16** | failed / **17** |

- **3D diode blocks at 0.5 V (`[.scale_carrier_blocks]`, Release, pinned to the P cores), iterations and the
  symmetrized AMG's setup + solve time:**

| Size (unknowns) | ψ (as is = symmetrized) | electrons as is | electrons symmetrized | holes as is | holes symmetrized |
|---|---|---|---|---|---|
| 20³ (24k) | 9 | failed | **10** (0.012 s) | failed | **11** (0.013 s) |
| 30³ (81k) | 11 | 750 | **13** (0.047 s) | failed | **14** (0.046 s) |
| 40³ (192k) | 11 | failed | **13** (0.13 s) | failed | **17** (0.15 s) |
| 50³ (375k) | 12 | failed | **15** (0.29 s) | failed | **19** (0.35 s) |
| 60³ (648k) | 13 | failed | **19** (0.64 s) | failed | **30** (0.90 s) |
| 70³ (1.03M) | 12 | failed | **21** (1.12 s) | failed | **35** (1.66 s) |

  ("failed": 1000 products without convergence.) The symmetrization turns blocks the AMG cannot solve into blocks it
  solves in 10–35 iterations. No unknown was clamped at 0.5 V (ln d within ±8.03).
- **MD1's go/no-go: no-go as written.** MG5's 1.5× on the blocks alone: the fixture meets it (electrons 11 → 13,
  1.18×; holes 15 → 17, 1.13×, from 48²), the 3D series does not: electrons 10 → 21 (**2.1×**), holes 11 → 35
  (**3.2×**), from 24k to 1.03M. The growth sits at 60³–70³, where the coarsest level grows (holes 766 → 2041,
  electrons 805 → 2376: the coarsening slows on the high-field rows, as the risks foresaw). Not loosened. Slice 1
  waits for the owner.
- For scale: at 70³ one symmetrized V-cycle solve of a whole carrier block (35 iterations) takes 1.7 s, against
  68.5 s for one GMRES-ILUT solve of the coupled system (Unit 58 slice 3); that is an observation, not a criterion.
- Mutants (MG9) not yet run for this slice; the full suites not run (owner's word).

**Owner ruling (2026-11-20, after slice 0): "2"**, investigate the growth first (slice 0b) before any slice 1.

**Slice 0b plan (fixed before any measurement):** the 3D blocks' coarsest levels lie above `coarse_size` (1000) with
fewer than `max_levels` levels, so the coarsening **stalled** (a level shrinking by less than 10%) and that level is
only smoothed (slice 2's rule), not solved; the ψ block never stalls (coarsest 268–883, factorized).

| # | Item |
|---|---|
| MB1 | **Diagnose (no product change):** per level of the symmetrized electron and hole blocks at 40³–70³: rows, aggregates, rows in no aggregate, rows with no strong neighbour, whether the last level stalled, and the share of its rows that are diagonally dominant by 10× or more (\|a_ii\| ≥ 10 Σ_j≠i \|a_ij\|). |
| MB2 | **Candidates (this list only):** (a) the strength threshold θ at 0.04, 0.02 and 0 (default 0.08); (b) a stalled coarsest level factorized by SparseLU when it has at most 20 × `coarse_size` unknowns (smoothed as now above that), a new `AmgSettings` flag, off by default; (c) (b) with (a)'s best θ. Each on the electron and hole blocks, 20³–70³, as `[.scale_carrier_blocks]`; the ψ block reported. |
| MB3 | **Selection:** a candidate passes if **both** carrier blocks grow at most 1.5× from 24k to 1.03M unknowns (MG5's target, on the blocks). Among those passing, the least 70³ setup + solve time of the two blocks; a tie (within 5%) goes to the fewest changes from the defaults. If none passes, the best is recorded and the owner asked again. |
| MB4 | **Unchanged elsewhere:** `amg` and `cpr` keep their settings and bits (MG1); a chosen candidate applies only to the carrier-block hierarchies of slice 1 (MD3) and is recorded there. |

**Slice 0b results (2026-11-20).** `AmgSettings::factorize_stalled` (off by default; `amg::stalled_factor` = 20) and
the hidden `[.carrier_0b]` in `solver_profile_test.cpp` (Release, pinned to the P cores; the symmetrized AMG as
slice 0, GMRES(30) to 1e-8). `nitcad_linalg_test` 108 cases pass with the flag (MB4).
- **MB1, the diagnosis (θ 0.08, smoothed; electrons, holes the same picture):** level 0 aggregates well (0.17 of
  its rows; the 9800 rows in no aggregate at 70³ are the contacts' identity rows). From level 2 on, most rows have
  **no strong neighbour** (70³ electrons: 2031 of 16 175 at level 2, 2146 of 3566 at level 3, 2151 of 2376 at level
  4), yet only 0.9–9.6% of the last level's rows are 10× diagonally dominant. So the coarse operators are not weakly
  coupled; their couplings are spread over many small Galerkin entries, none of which reaches θ √(a_ii a_jj). The
  coarsening ratio rises to 0.65–0.93 and stalls (50³–70³ electrons, 60³–70³ holes).
- **MB2, iterations (the setup + solve time in seconds at 70³):**

| Block, candidate | 20³ | 30³ | 40³ | 50³ | 60³ | 70³ | 70³ / 20³ | worst / 20³ |
|---|---|---|---|---|---|---|---|---|
| electrons, θ 0.08 (as slice 0) | 10 | 13 | 13 | 15 | 19 | 21 (1.26) | 2.10 | 2.10 |
| electrons, θ 0.04 | 12 | 13 | 13 | 20 | 24 | 16 (1.10) | **1.33** | 2.00 |
| electrons, θ 0.02 | 13 | 15 | 17 | 17 | 28 | 24 (1.39) | 1.85 | 2.15 |
| electrons, θ 0 | 13 | 15 | 14 | 16 | 19 | 23 (1.36) | 1.77 | 1.77 |
| electrons, θ 0.08, factorized | 10 | 13 | 13 | 14 | 18 | 19 (1.32) | 1.90 | 1.90 |
| electrons, θ 0.04, factorized | 12 | 13 | 13 | 20 | 24 | 16 (1.15) | **1.33** | 2.00 |
| holes, θ 0.08 (as slice 0) | 11 | 14 | 17 | 19 | 29 | 35 (1.93) | 3.18 | 3.18 |
| holes, θ 0.04 | 13 | 14 | 16 | 17 | 24 | 16 (1.08) | **1.23** | 1.85 |
| holes, θ 0.02 | 15 | 15 | 17 | 18 | 27 | 25 (1.48) | 1.67 | 1.80 |
| holes, θ 0 | 13 | 15 | 14 | 18 | 19 | 17 (1.13) | 1.31 | 1.46 |
| holes, θ 0.08, factorized | 11 | 14 | 17 | 19 | 19 | 22 (1.53) | 2.00 | 2.00 |
| holes, θ 0.04, factorized | 13 | 14 | 16 | 17 | 24 | 16 (1.25) | **1.23** | 1.85 |

  (θ 0.02 and 0 factorized equal their smoothed rows within one iteration: at those θ nothing stalls.) The ψ block:
  9–15 iterations under every candidate, 70³ 12–14.
- **MB3, the selection.** Read as written (24k against 1.03M), two candidates pass on both blocks: **θ 0.04**
  (electrons 1.33×, holes 1.23×) and θ 0.04 factorized (the same counts; (c), since 0.04 is (a)'s best). θ 0.04
  smoothed takes 2.18 s at 70³ for the two blocks against 2.40 s factorized, so by MB3 the choice is **θ = 0.04
  for the carrier-block hierarchies, stalled levels smoothed** (one change from the defaults).
- **But the growth is not monotone:** at 60³ θ 0.04 needs 24 iterations on both blocks, 2.0× and 1.85× the 20³
  counts, and 70³ falls back to 16. No candidate keeps every size within 1.5× (the best, holes at θ 0, peaks at
  1.46× but electrons at θ 0 reach 1.77×). MB3 compares the end points only, so θ 0.04 passes as written; whether
  the 60³ peak counts against MG5 is put to the owner, not decided here.
- **What the diagnosis points to (not in MB2's list, not run):** a strength threshold that falls with the level
  (θ_l = θ · ½^l, attributed from memory to Vaněk–Mandel–Brezina 1996, not verified: their paper's text was not reachable), since the stall comes from the coarse levels' spread
  couplings, not from level 0.

**Owner ruling (2026-11-20, after slice 0b): "Search in research papers", then "yes"** to a slice 0c with three
candidates from the literature and a stricter criterion. Sources read (abstracts only, full texts not reached): Lin,
Shadid, Tuminaro, Sala, IJNME 2010 (Petrov–Galerkin smoothed aggregation on drift-diffusion; not a candidate here,
it replaces MD2); Olson, Schroder, Tuminaro, NLAA 17 (2010) 713–733 (the evolution strength measure); sparsified /
non-Galerkin smoothed aggregation (Galerkin fill spreads coarse couplings); Clees and Stüben, SAMG device
simulation (shared interpolation over the unknowns).

**Slice 0c plan (fixed before any measurement):**

| # | Item |
|---|---|
| MC1 | **Candidates (this list only), each a new `AmgSettings` option, off by default (`amg` and `cpr` unchanged, MB4):** (1) a **falling threshold**: level l uses θ · ½^l, with θ 0.08 and 0.04; (2) the **evolution measure, simplified** (the paper's text was not reached, so this is our reading of the abstract, not its algorithm): z = (I − ω D⁻¹A)^4 δ_i with ω = 1 / max_r Σ_j \|a_rj / a_rr\| (Gershgorin), m_ij = \|1 − (z_j b_i) / (z_i b_j)\| with b the level's near-null vector; j strong for i if m_ij ≤ ε · min_k m_ik (over i's off-diagonal neighbours, the minimum floored at 1e-12), with ε 2 and 4; (3) a **row-relative threshold** (classical AMG's test, which keeps a row's largest couplings however spread): j strong for i if \|a_ij\| ≥ θ max_{k≠i} \|a_ik\|, with θ 0.25 and 0.5. Strength stays symmetrized and pass 2 still ranks by \|a_ij\| / √(a_ii a_jj). Stalled levels smoothed (slice 2's rule). |
| MC2 | **Measured** as slice 0b on the symmetrized electron and hole blocks (and the ψ block, reported), 20³–70³, plus the slice 0b default (θ 0.08) and θ 0.04 again as references. |
| MC3 | **Criterion (stricter than MB3, never loosened):** a candidate passes if, on **both** carrier blocks, **every** size's iterations are at most 1.5× that block's 20³ count. Among those passing, the least 70³ setup + solve time of the two blocks; within 5%, the simplest (in the order 1, 3, 2). If none passes, recorded, and slice 1 is put to the owner with the best result. |
| MC4 | Setup time reported for every candidate (the evolution measure's cost grows with its four steps' support). |

**Slice 0c results (2026-11-20).** `StrengthMeasure` (`scaled`, `row_max`, `evolution`) and `strength_decay` in
`AmgSettings` (defaults `scaled`, 1: `amg` and `cpr` unchanged; `nitcad_linalg_test` 108 cases pass); hidden
`[.carrier_0c]` (Release, pinned to the P cores). Iterations (setup + solve in seconds at 70³):

| Block, candidate | 20³ | 30³ | 40³ | 50³ | 60³ | 70³ | worst / 20³ |
|---|---|---|---|---|---|---|---|
| electrons, θ 0.08 (reference) | 10 | 13 | 13 | 15 | 19 | 21 (1.14) | 2.10 |
| electrons, θ 0.04 (reference) | 12 | 13 | 13 | 20 | 24 | 16 (1.01) | 2.00 |
| electrons, falling 0.08 | 10 | 10 | 16 | 18 | 19 | 14 (0.83) | 1.90 |
| electrons, falling 0.04 | 12 | 14 | 16 | 18 | 29 | 21 (1.23) | 2.42 |
| electrons, evolution ε 2 | 16 | 19 | 27 | 28 | 35 | 37 (8.98) | 2.31 |
| electrons, evolution ε 4 | 15 | 19 | 23 | 28 | 31 | 30 (8.69) | 2.07 |
| electrons, row-relative 0.25 | 11 | 12 | 11 | 21 | 18 | 18 (1.09) | 1.91 |
| electrons, row-relative 0.5 | 10 | 13 | 15 | 17 | 25 | 25 (1.28) | 2.50 |
| holes, θ 0.08 (reference) | 11 | 14 | 17 | 19 | 30 | 35 (1.65) | 3.18 |
| holes, θ 0.04 (reference) | 13 | 14 | 16 | 17 | 24 | 16 (0.95) | 1.85 |
| holes, falling 0.08 | 11 | 12 | 17 | 21 | 21 | 18 (0.96) | 1.91 |
| holes, falling 0.04 | 13 | 14 | 17 | 18 | 27 | 19 (1.06) | 2.08 |
| holes, evolution ε 2 | 13 | 16 | 14 | 18 | 24 | 22 (5.47) | 1.85 |
| holes, evolution ε 4 | 13 | 16 | 14 | 17 | 21 | 27 (5.64) | 2.08 |
| holes, row-relative 0.25 | 11 | 12 | 11 | 21 | 18 | 13 (0.88) | 1.91 |
| holes, row-relative 0.5 | 11 | 15 | 16 | 18 | 23 | 22 (1.18) | 2.00 |

  The ψ block: 9–16 iterations under every candidate (falling 0.08: 9–10 at every size, its best).
- **MC3: no candidate passes.** Every one exceeds 1.5× at some size on at least one block. The best worst-case
  ratios are falling 0.08 and row-relative 0.25 (1.90–1.91); of these, **falling 0.08** is the fastest at 70³
  (1.79 s for the two blocks against 1.97 s) and the simplest, so it is **the best result**. Its end points are
  1.40× (electrons 10 → 14) and 1.64× (holes 11 → 18); its peaks are at 40³–60³ (16–21).
- **The counts are not monotone in size** under any test (peaks at 50³ or 60³ fall back at 70³), so the growth is
  not a smooth function of h: the aggregates' alignment with the junction and the high-field rows changes from mesh
  to mesh. Changing the strength test moves the peak; it does not remove it.
- **The simplified evolution measure is the worst** on electrons and costs 5–9 s of setup at 70³ (its 4-step
  supports), 5–8× the others. Our reading of the abstract is not the paper's algorithm; this rules out only the
  reading, not the method.
- **Per MC3, slice 1 is put to the owner** with the best result (falling 0.08). Not loosened.

**Owner ruling (2026-11-20, after slice 0c): "start slice 1"**, the recommended option: build `cpr_amg` with the
falling threshold (the `amg_*` θ halved per level, `strength_decay` ½) on every stage-2 field hierarchy, judged on
MG2, MG3, MG6–MG8; MG5's miss on the blocks alone (worst 1.90–1.91×) stays recorded. Stage 1's AMG on the reduced
potential system keeps `cpr`'s settings.

**Slice 1 (MD3, MD4), 2026-11-20: `preconditioner = cpr_amg`, built.** `CprBackend` gains `Stage2::block_amg`:
stage 1 as `cpr`, then stage 2 one forward block Gauss–Seidel sweep of a `SymmetrizedAmg` (falling θ, `strength_decay`
½) per carrier field on r2 = v − A z1, then a Gauss–Seidel step for each unknown in no node. Wired through
`LinearSolver` (real only, needs labels and `backend = gmres`), `newton_solve`'s labels, trace and weighted current
(the Laplace problem takes `amg`), the deck (`preconditioner = cpr_amg`) and the app's vocabulary. Release only (no
Debug by standing rule); targeted runs: linalg (112 cases), deck, app, solve `[gmres]` and the retry case, all passed.

Two findings changed MD3 as proposed, both recommended options taken under the standing rule:

- **The potential block is not swept (Sweep::carriers).** MD3 swept every field, ψ first by AMG on A_ψψ alone. On
  the drift-diffusion fixture that diverged from 48² up (GMRES 1000 iterations, every bias) while sweeping only the
  carrier blocks converged everywhere: A_ψψ lacks the screening the carriers add, and stage 1 already owns ψ. Probe
  `[.cpr_amg_sweep]`, GMRES iterations (cpr_amg carriers only / MD3 all fields / cpr):

  | fixture | bias 0 | bias 20 |
  | --- | --- | --- |
  | 24² | 2 / 4 / 30 | 37 / 38 / 28 |
  | 48² | 18 / failed / 101 | 38 / failed / 54 |
  | 96² | 22 / failed / failed | 43 / failed / 300 |
  | 192² | 21 / failed / failed | 41 / failed / failed |

  (A further variant, stage 1 again on the sweep's residual, gave the carriers-only counts exactly; not kept.)
- **A Jacobian the sweep leaves GMRES short on is retried once with ILUT as stage 2 (`cpr`), strict**, counted in
  `BackendCounts::stage2_retries`; that factorization keeps ILUT until the next. Found by MG3's small-signal
  junction case (p+n 1e19/1e16, −2 V then 0.6 V): the first iterate of the 0.6 V point, with the new bias stamped
  on the contacts, failed in 1000 iterations where `cpr` takes 5. There even a direct SparseLU of the hole block
  alone reached only backward error 1.25e-1 (1.5e-16 before the stamp); row equilibration fixed the bare LU
  (4.7e-15) but not the symmetrized solve (14 iterations for an exact one-level solve) and doubled the fixture's
  counts (27,648 unknowns at bias 20: 43 → 87), so it was not kept. Two guards on the scaling were tried and
  dropped: the overall span of ln d (56 at the solved −2 V too, which needs no help) and the largest jump across
  a usable pair (2.9 at both states). The cause is not pinned down; the retry is measured, not explained.
  Test: the stamped step needs one retry, the solved point none.

Results against the gates so far: **MG2 met** (12,288 unknowns, 61 iterations, the same bits on a repeat and at 2,
4, 6, 8 threads). **MG3 met** on the cases `cpr` has (diode J(0.5 V) to 1e-10 and the forward sweep, worst 0.000;
MOS C-V; transient, worst 7.0e-21 of the largest current; small-signal junction 0.133 of the bound, with the retry).
The fixture against `cpr` (GMRES iterations, cpr_amg / cpr): 1,728 unknowns 2 / 24 and 36 / 14; 6,912: 27 / 60 and
39 / 36; 27,648: 40 / 330 and 43 / 120 (bias 0 and 20). MG6–MG8 (3D series, FinFET against PARDISO) are MD5's
measurements, next.

**MG9 (mutants), 2026-11-20: 11 of 12 caught** (script `C:\Users\disha\ssrun\mut63.sh`, the later ones run one at a
time, Release, pinned to 4 cores by the owner's word). Caught: M2 every field swept, M3 no retry, M4 the retry keeping
the sweep, M5 strength decay 1, M6 stage 2 ILUT, M7 the deck refusing `cpr_amg`, M8 the sweep's result scattered at
half weight, M9 labels not required, M10 Newton giving no labels, M11 the unknowns in no node given zero, M12 complex
accepted. Two tests were too weak and were strengthened:

- **M11 survived because the retry hid it.** With the bordering unknown left at zero the block sweep misses and the
  ILUT retry converges (13 iterations against 19), so a solve that never used the sweep passed. The `cpr_amg` tests
  now require no stage-2 retry (`Solved::stage2_retries`).
- **M8 survived because only the answer was checked.** The fixture test now bounds GMRES iterations at the counts
  above plus about a tenth (4, 40; 30, 43; 44, 47); M8 needs 30, 49; 45, 56.

**M1 (the sweep without the field-to-field coupling, block Jacobi) survives, and is a finding, not a test gap.** The
fixture's carrier blocks have no electron–hole entries, so M1 is equivalent there. With a recombination term r on
each interior node (r on both diagonals and both off-diagonals, `with_recombination`; probe
`[.cpr_amg_recombination]`), GMRES iterations, forward sweep (as built) / M1 Jacobi / `cpr` (– = missed in 1000):

| unknowns, bias | r = 0.5 | r = 0.05 | r = 0.005 |
|---|---|---|---|
| 1,728, 0 | 90 / 77 / – | 21 / 25 / – | 13 / 10 / 109 |
| 1,728, 20 | 142 / 75 / 780 | 50 / 48 / 209 | 41 / 40 / 57 |
| 6,912, 0 | 266 / 180 / – | 56 / 59 / – | 33 / 28 / – |
| 6,912, 20 | 259 / 120 / – | 58 / 53 / – | 40 / 40 / – |
| 27,648, 0 | 900 / 660 / – | 300 / 120 / – | 55 / 59 / – |
| 27,648, 20 | 720 / 531 / – | 210 / 148 / – | 49 / 52 / – |

`cpr_amg` solves every coupled case and `cpr` almost none, but both orderings grow with r and with size (55 → 900
at 27,648), and the coupling the forward sweep carries makes it worse, not better, from r = 0.05 up. The owner
chose (2026-11-20, "choose recommended") a coupled carrier stage, on this table.

**The table's coupling is not physical, and SRH's shape reverses the finding (2026-11-20, before any coupled stage
was built).** A constant r on both rows is not what recombination stamps: SRH, R = (np − 1) / (n + p + 2) (n_i = n_1
= p_1 = 1, unit lifetimes) times r (`with_srh`), puts r ∂R/∂n on the electron diagonal and the hole row's electron
column and r ∂R/∂p on the hole diagonal and the electron row's hole column, so each row's coupling follows the other
carrier's density, as the symmetrization's scaling does. Same probe, GMRES iterations, forward sweep (as built) /
M1 Jacobi / `cpr`:

| unknowns, bias | r = 1000 | r = 100 | r = 10 | r = 1 | r = 0.1 |
|---|---|---|---|---|---|
| 1,728, 0 | 5 / 9 / 22 | 4 / 6 / 22 | 4 / 7 / 22 | 4 / 5 / 17 | 3 / 4 / 23 |
| 1,728, 20 | 162 / 271 / 11 | 58 / 90 / 13 | 35 / 37 / 14 | 39 / 37 / 14 | 37 / 36 / 15 |
| 6,912, 0 | 30 / 34 / 53 | 28 / 28 / 61 | 27 / 28 / 55 | 27 / 27 / 58 | 27 / 27 / 60 |
| 6,912, 20 | 270 / 540 / 36 | 97 / 171 / 34 | 41 / 54 / 31 | 39 / 38 / 34 | 39 / 39 / 36 |
| 27,648, 0 | 46 / 44 / – | 45 / 43 / – | 40 / 40 / – | 41 / 37 / – | 38 / 39 / – |
| 27,648, 20 | 660 / 420 / 420 | 180 / 326 / 300 | 55 / 100 / 267 | 49 / 49 / 390 | 44 / 41 / 321 |

Up to r = 10 the counts are the uncoupled fixture's (40–55 at 27,648), and the forward sweep is as good as Jacobi
or better (55 against 100 at 27,648, bias 20, r = 10): the coupling it carries helps. Only strong recombination
under bias (r ≥ 100, where h exceeds the diffusion length by a factor of 10) grows, by up to 15 times, and `cpr`
grows there as well. M1 stays a survivor: Jacobi is the weaker form, but no count bound separates the two
reliably (they cross at 27,648, bias 20, r = 1000). The owner then took the recommended option (2026-11-20,
"choose recommended"): keep the forward sweep, measure MD5, and build a coupled carrier stage only if MD5's devices
show the strong regime.

**MD5, 2026-11-20: MG6 and MG8 met, MG7 met on the FinFET sweep but missed on the diode series' threads, MG5
missed.** Release, pinned to the six P cores (mask 3C3; the owner's four-core cap of earlier the same day was
lifted before these runs). The 3D diode series is `[.scale_cpr_amg]` (log `C:\Users\disha\ssrun\md5series.log`):
one Newton step's linear solve, the Jacobian at 0.5 V on the equilibrium potential, a fresh solver, strict.
Seconds (build + solve), GMRES products and held memory:

| size, unknowns | cpr_amg, 1 thread | cpr_amg, 6 threads | GMRES-ILUT | GMRES-CPR |
|---|---|---|---|---|
| 20³, 24,000 | 0.06 + 0.06, 22, 34 MB | 0.06 + 0.06 | 0.56 + 0.16, 40, 46 MB | 0.57 + 0.22, 42, 57 MB |
| 30³, 81,000 | 0.23 + 0.37, 28, 122 MB | 0.22 + 1.61 | 2.97 + 1.22, 60, 156 MB | 3.21 + 1.36, 60, 192 MB |
| 40³, 192,000 | 0.60 + 1.47, 42, 269 MB | 0.58 + 1.49 | 7.82 + 3.17, 90, 379 MB | 6.98 + 2.98, 71, 471 MB |
| 50³, 375,000 | 1.00 + 3.31, 52, 526 MB | 1.19 + 1.52 | 15.23 + 8.61, 108, 735 MB | 16.13 + 9.15, 96, 905 MB |
| 60³, 648,000 | 1.81 + 7.25, 71, 918 MB | 1.97 + 3.69 | 29.03 + 21.21, 116, 1,270 MB | 29.70 + 21.43, 112, 1,570 MB |
| 70³, 1,029,000 | 3.71 + 17.83, 83, 1,437 MB | 3.50 + 7.37 | 46.76 + 35.36, 142, 2,017 MB | 48.37 + 38.81, 138, 2,487 MB |

- **MG6 met:** at one thread the linear time is 0.19× GMRES-ILUT's at 40³ and 0.26× at 70³ (target ≤ 0.5×), and
  0.15–0.26× at every size; GMRES-CPR is no faster than ILUT at any size.
- **MG8 met:** held memory 0.71× GMRES-ILUT's at 40³ and at 70³ (269 against 379 MB, 1,437 against 2,017).
- **MG7, the series: missed.** Six threads against one: 1.0× at 40³, 1.98× at 70³ (the solve alone 2.42×), against
  ≥ 2.5×; the build (aggregation, symmetrization) is serial. At 30³ six threads are slower (solve 1.61 s against
  0.37 s, as `[.scale_amg]`'s 70³ six-thread anomaly in Unit 58); not explained. The bits agree (MG2): the same
  GMRES products at 1 and 6 threads at every size.
- **MG7, the FinFET sweep: met.** `decks/finfet_3d.deck` (45,375 nodes, 17 points), copies with `[solver]` keys
  only (`C:\Users\disha\ssrun\md5`), both at 6 threads: PARDISO 526 s (17:16:47–17:25:33), `gmres` + `cpr_amg`
  224 s (17:25:33–17:29:17), 2.35× faster. The drain currents agree at every point to 7.6e-11 relative (worst
  absolute 3.7e-20 A, at Vg = 0), within any terminal-current resolution. `cpr_amg` took more Newton iterations
  at Vg = 0.30–0.40 V (point 6 retried: 32, then 38; points 7 and 8 20 and 18, against PARDISO's 7 and 8), and
  PARDISO retried point 1 (33, then 41) where `cpr_amg` took 11: the strong-recombination regime of the probe did
  not show (no point failed and no linear solve fell back), so no coupled carrier stage is built.
- **MG5 missed:** GMRES products grow 3.8× from 24,000 to 1,029,000 unknowns (22 → 83; target ≤ 1.5×, `cpr` 3.5×
  here) and on the fixture from 24² to 96² 1.19× at bias 20 (36 → 43) but from 2 to 40 at bias 0. Slices 0b–0c
  had already found no strength test that meets 1.5× at every size.
- A first run at four threads (the cap then in force; `md5series_4threads.log`) gave the same one-thread figures
  within timing noise but 77 GMRES products at 60³ against 71 here: the series' equilibrium comes from PARDISO at
  6 threads, which is not bit-reproducible between runs, so the Jacobian differs in its last bits. Within each run
  the thread counts agree.

**Close-out, 2026-11-20 (owner: "choose recommended", close with the misses recorded as limitations):**
- **ARCHITECTURE.md:** 6.10 gains "As built (Unit 63)" (the symmetrization, `cpr_amg`, the measurements and the
  limits), and the 11.3 row says what Unit 63 built; its next steps are the iteration growth, thread scaling (the
  serial AMG build) and a coupled electron–hole stage for strong recombination.
- **Unit 63's results:** MG2 met (bits at 1 to 8 threads); MG3 met; MG4 met (slice 0); MG6 met; MG8 met; MG9 met
  (11 of 12 caught, M1 recorded); MG10 met (labels stay plain integers, no Eigen in public headers, no new
  dependency); MG7 met on the FinFET sweep, missed on the diode series' threads; MG5 missed (recorded, not
  loosened). **MG1 met:** against `origin/main` 2065afc (built in a temporary worktree, Release), every printed count
  line of the `[gmres]`/`[amg]` cases and of the whole linalg binary (25 lines) appears unchanged on the branch,
  which only adds the `cpr_amg` lines. The full Release and Debug suites have not run (owner's word needed; no Debug build by
  standing rule). Nothing pushed.
### Unit 66 — Philips unified (Klaassen) mobility, opt-in (proposed 2026-11-23; on `physics/klaassen-mobility`)

**Owner request (2026-11-23):** "Implement Philips unified (Klaassen) model in my NiTCAD", then "choose recommended and
start". Klaassen's model was one of the top gaps of the 2026-11-16 review (Shi 2026). It gives one description of the
majority and minority mobilities of both carriers in silicon: lattice scattering, donor and acceptor scattering with
screening by the carriers, electron–hole scattering and dopant clustering at high doping. It is the bipolar, power
and solar-cell default of the commercial tools; Caughey-Thomas and Masetti give a minority carrier the majority
mobility of the same total doping, which is wrong by up to a factor of about 2 in heavily doped emitters.

**Standing conditions.** Opt-in. Unselected, everything is bit-identical to `main`: residuals, Jacobian patterns and
values, run identities, result files, the Unit 40/59 traces and decks. Criteria are fixed before code and never
loosened after a run (TVM2). Release and Debug suites with `/W4 /WX`, run only on the owner's word. Commits stay
local; no push, PR or merge without the owner's word. The approach (staged, KQ1) and KQ2–KQ7 were ruled as
recommended under the owner's standing instruction (2026-11-11) and the owner's "choose recommended and start".

**Sources (read 2026-11-23; the 1992 papers themselves were not read, no copy in hand).** D. B. M. Klaassen,
Solid-State Electron. 35, 953 (1992) and 961 (1992), as reproduced by: (S1) Sandia Charon User Manual v2.2,
SAND2022-7653, section 7.4, eqs (62)–(74) and Tables 7-8, 7-9; (S2) DEVSIM, "Mobility modeling" (2009), its model
code; (S3) COMSOL Semiconductor Module User's Guide 6.3/6.4, "Klaassen Unified Mobility Model", parameter list;
(S4) Silvaco ATLAS User's Manual (an old edition), eqs 3-158–3-183. Copies stay in `C:\Users\disha\papers\klaassen\`
(not committed). Where they disagree:
- r1, r5: S1 and S3 give 0.7643 and −0.01552; S4 gives 0.7463 and −0.8552. Taken: S1/S3 (two independent sources;
  S4's r1 is a digit transposition).
- s1: S1 prints 0.89223, S3 and S4 0.89233. Taken: 0.89233.
- Electron set: S4's defaults (1417, 52.2, 9.68e16, 0.68) are Klaassen's arsenic set (S1 Table 7-8); S1 and S3
  give phosphorus (1414, 68.5, 9.20e16, 0.711) as the electron default. Taken: phosphorus (KQ4).
- P_BH: S1, S2 use n + p; S4 uses n alone for electrons and p alone for holes. Taken: n + p (S1, S2).
- Clustering: S1 eq (74) and S2 apply it to the dopant densities; S4 only inside P_CW. Taken: S1 (KQ5).
- P_CW: S2 writes 3.97e13 (1/(Z³ N_sc)) ^(2/3) (T/300)² with Z the clustering factor; S1 3.97e13 N_sc^(−2/3) with
  the temperature factor not legible; S4 3.97e13 (1/(Z³ N_D))^(2/3) (T/300)². Taken: 3.97e13 N_sc^(−2/3) (T/300)²,
  Z = 1 (singly ionised dopants; the clustering already enters through N*), KQ5.
- G(P) minimum: S1 holds G at its minimum for P below the minimiser (G otherwise turns up and can go negative for
  small P). Taken (KQ6).

**Model.** For carrier i (e or h), the other carrier j, lattice temperature T [K], N_D, N_A [cm⁻³], n, p [cm⁻³]:

  1/μ_i = 1/μ_i,L + 1/μ_i,DAP                                                        (Matthiessen)
  μ_i,L = μ_max,i (T/300)^(−θ_i)
  μ_i,DAP = μ_i,N (N_i,sc / N_i,sc,eff) (N_ref,i / N_i,sc)^α_i + μ_i,c (n + p) / N_i,sc,eff
  μ_i,N = μ_max,i² / (μ_max,i − μ_min,i) · (T/300)^(3α_i − 1.5)
  μ_i,c = μ_max,i μ_min,i / (μ_max,i − μ_min,i) · (300/T)^0.5
  N_D* = N_D (1 + 1 / (c_D + (N_ref,D / N_D)²)),  N_A* likewise with c_A, N_ref,A (N* = 0 at N = 0)
  N_e,sc = N_D* + N_A* + p,             N_h,sc = N_A* + N_D* + n
  N_e,sc,eff = N_D* + G(P_e) N_A* + p / F(P_e),   N_h,sc,eff = N_A* + G(P_h) N_D* + n / F(P_h)
  P_i = [f_CW / P_CW,i + f_BH / P_BH,i]⁻¹,  P_CW,i = 3.97e13 N_i,sc^(−2/3) (T/300)²,
                                             P_BH,i = 1.36e20 / (n + p) · (m_i/m0) (T/300)²
  G(P) = 1 − s1 / [s2 + ((m0/m_i)(T/300))^s4 P]^s3 + s5 / [((m_i/m0)(300/T))^s7 P]^s6, held at its minimum below
         the minimiser P_min,i (found once per (i, T) by Brent on the monotone derivative)
  F(P_i) = (r1 P^r6 + r2 + r3 m_i/m_j) / (P^r6 + r4 + r5 m_i/m_j)

Evaluated as 1/μ_i,DAP = N_i,sc,eff / [μ_i,N N_ref,i^α N_i,sc^(1−α) + μ_i,c (n + p)], which is 0 (μ = μ_i,L) when
N_sc = 0 and n + p = 0; P = ∞ gives G = 1, F = r1 by their limits.

Parameters (silicon): electrons (phosphorus) μ_max 1414, μ_min 68.5 cm²/(V s), N_ref 9.20e16 cm⁻³, α 0.711,
θ 2.285; holes (boron) 470.5, 44.9, 2.23e17, 0.719, θ 2.247; c_D 0.21, c_A 0.50, N_ref,D 4.0e20, N_ref,A 7.2e20;
m_e/m0 = 1, m_h/m0 = 1.258; f_CW 2.459, f_BH 3.828; s1–s7 0.89233, 0.41372, 0.19778, 0.28227, 0.005978, 1.80618,
0.72169; r1–r6 0.7643, 2.2999, 6.5502, 2.3670, −0.01552, 0.6478 (S1, S3).

**Design (KD1–KD9).**

| ID | Decision | Reason | Limitation |
|---|---|---|---|
| KD1 | Physics: `physics/klaassen_mobility.hpp/.cpp`, `KlaassenParameters` (the set above, `klaassen_silicon`) and `klaassen_mobility(p, carrier, N_D, N_A, n, p, T)` returning μ and ∂μ/∂n, ∂μ/∂p, ∂μ/∂T, all analytic. Preconditions as Masetti's. | One function serves both variants and every assembly path; analytic partials for Newton (ARCHITECTURE section 5). | — |
| KD2 | Two deck values per carrier (KQ1, KQ2): `electron_mobility` / `hole_mobility` = `klaassen` (n, p the solution's) or `klaassen_equilibrium` (n, p the neutral equilibrium of the node's doping: n₀ = N/2 + √((N/2)² + n_i²), p₀ = n_i²/n₀ with N = N_D − N_A, n_i the material's at T; static, like Masetti). Both carriers may differ. | `klaassen_equilibrium` is the drop-in (every path, no Jacobian change); `klaassen` is the model's purpose (high injection). | — |
| KD3 | Silicon only: every semiconductor region must carry silicon's Caughey-Thomas mobility set, as Masetti (Units 34, 38); refused (`invalid_input`) otherwise. `doping_mobility off` makes N_D = N_A = 0 (μ then from carriers and lattice only). | The fits are silicon's. | Other materials: follow-up. |
| KD4 | `klaassen_equilibrium` (slice 0): the node mobilities of drift_diffusion.cpp's precompute, electrothermal's node table (μ and ∂μ/∂T, as `masetti_mobility_slope`; n₀, p₀ re-evaluated at the node's T), energy balance and impact ionization read it unchanged. | No pattern change; bit-identical paths for everything else. | Carrier–carrier scattering at its equilibrium value only. |
| KD5 | `klaassen` (slice 1): the node mobilities become functions of the node's n, p. Each Newton evaluation recomputes μ_n, μ_p and their n, p partials per semiconductor node; the edge low-field mobility stays the harmonic mean of its ends, so its partials reach n_a, p_a, n_b, p_b. The isothermal edge Jacobian gains the cross entries (electron row – p columns, hole row – n columns, both ends), added only when `klaassen` is selected. Field (Canali) and surface mobility chain through d μ / d μ₀ as today. | Exact Newton (quadratic convergence) is the project's rule; a lagged mobility would need outer iterations and change convergence behaviour. | Pattern grows 4 entries per edge direction when on. |
| KD6 | Refused with `klaassen` (KQ3; `invalid_input`, before solving): electrothermal, energy balance, density gradient, impact ionization, quantum `mlda`. `klaassen_equilibrium` is allowed with all of them. | Each of these keeps its own node-mobility table or derivative chain; carrier-dependent versions are follow-ups. | Recorded follow-ups. |
| KD7 | Run records: `models.electron_mobility` gains 2 (`klaassen`) and 3 (`klaassen_equilibrium`); `models.hole_mobility` gains 3 and 4 (2 is `automatic`). Written only when not the default, as now. | Run identities of existing decks unchanged. | — |
| KD8 | Deck vocabulary, app field list and docs: the two values added to `deck_vocabulary.cpp` with their help text; ARCHITECTURE.md section 6 as-built block and the 11.3 row. | The app's deck editor offers the values. | — |
| KD9 | Example: none new; the bundled diode deck gets a `klaassen` variant in a test only (KC7). | No new bundled file without the owner (the tree is the owner's). | — |

**Questions (ruled as recommended, 2026-11-23).**

| ID | Question | Options | Ruling |
|---|---|---|---|
| KQ1 | Scope | (a) equilibrium-carrier (doping-only) Klaassen only; (b) full carrier-dependent everywhere at once; (c) staged: (a) in every path in slice 0, full carrier dependence in isothermal drift-diffusion in slice 1, the rest refused and recorded | **(c)** |
| KQ2 | How the variants are selected | (a) two deck values `klaassen`, `klaassen_equilibrium`; (b) one value and a separate switch | **(a)** |
| KQ3 | `klaassen` with electrothermal, energy balance, density gradient, impact ionization, MLDA | (a) refused, follow-ups; (b) silently fall back to equilibrium carriers | **(a)**: no silent model change |
| KQ4 | Electron parameter set | (a) phosphorus (S1, S3 default); (b) arsenic (S4 default); (c) per-region dopant choice | **(a)**, arsenic a recorded follow-up |
| KQ5 | Clustering and P_CW | (a) clustering on the dopants (S1 eq 74), P_CW with N_sc and Z = 1; (b) S2's Z³ inside P_CW as well; (c) S4's form | **(a)** |
| KQ6 | G(P) below its minimiser | (a) hold G at its minimum (S1); (b) the raw formula | **(a)** |
| KQ7 | Published curves | (a) no digitised paper data (papers not in hand); KC2–KC5 check limits, Masetti agreement for majority carriers and the minority > majority property; (b) the owner supplies the papers first | **(a)**, recorded as open |

**Criteria (KC1–KC10), fixed before code.**

| ID | Criterion |
|---|---|
| KC1 | Off: every existing test passes unchanged; the Unit 59 NG1 trace set byte-identical to `main`'s; run identities and Jacobian patterns unchanged. |
| KC2 | Limits: at N_D = N_A = n = p = 0, μ = μ_L to 1e-14 relative at 200, 300, 400 K; at 300 K μ_e = 1414, μ_h = 470.5. |
| KC3 | Majority carriers at 300 K, n = N_D (N_A = 0, p = n_i²/n) and p = N_A: within 10 % of Masetti's phosphorus (electrons) and boron (holes) fits (physics/mobility.hpp) at every decade 1e14–1e20 cm⁻³; the largest deviation recorded. |
| KC4 | Minority over majority at 300 K: μ_e in p-type (N_A, p = N_A) > μ_e in n-type (N_D = N_A) at 1e19 and 1e20 cm⁻³; likewise μ_h; the ratios recorded. |
| KC5 | Partials: ∂μ/∂n, ∂μ/∂p, ∂μ/∂T against central differences to 1e-6 relative over a grid (N_D, N_A in {0, 1e15, 1e18, 1e20}, n, p in {0, 1e10, 1e16, 1e19}, T in {250, 300, 400}); G continuous at P_min. |
| KC6 | `klaassen_equilibrium`: a PN diode IV (bundled 1D diode deck) converges at every bias; the Jacobian of a two-carrier state matches central differences of the residual to 1e-6 relative (existing FD-check helper), isothermal and electrothermal. |
| KC7 | `klaassen`: the same Jacobian FD check (cross entries included) to 1e-6; the diode IV converges with Newton steps per bias no more than the Masetti run's plus 2; in high injection (forward bias with n, p > 10 N_D in the base) the current differs from `klaassen_equilibrium`'s (carrier–carrier scattering acts), the difference recorded. |
| KC8 | Refusals: each KD3 and KD6 refusal names the region or model; the deck errors for unknown values list all values. |
| KC9 | Run records and restart: the new values recorded; a resumed `klaassen` transient continues bit for bit (Unit 61). |
| KC10 | Release and Debug suites (owner's word), NG1 traces vs `main`; `/W4 /WX` clean. |

**Slices.** 0: physics + `klaassen_equilibrium` in every path + deck/records (KC2–KC6, KC8 partly). 1: `klaassen`
carrier-dependent assembly and refusals (KC7–KC9). 2: ARCHITECTURE as-built + suites on the owner's word (KC1, KC10).

**Slice 0 record, physics (2026-11-23).** `physics/klaassen_mobility.hpp/.cpp` (forward-mode partials in n, p, T;
G's minimiser by 60 bisections of G' on ln P in [1e-8, 1e8]), `tests/physics/klaassen_test.cpp`. The formulas match
an independent transcription (golden-section minimiser) to 1e-10 at 3 × 8 × 8 × 5 × 5 points per carrier. KC2 met.
KC3 met: largest deviation from Masetti electrons 0.56 % (1e20), holes 4.64 % (1e18). KC4 met: minority/majority
at 1e19 electrons 1.640, holes 1.790; at 1e20 electrons 2.261, holes 2.542. KC5: two readings, stated here and
put to the owner as KQ8:
- Written into the test before its first run: the points with N_D = N_A = 0 and a zero carrier density are left
  out (not differentiable, the header says so).
- **Changed after the first run (TVM2: a loosening, not self-ruled).** The first run failed 675 of 13410 checks,
  all partials whose whole effect x ∂μ/∂x is under the central difference's rounding noise (about 5e-12 μ at
  h = 1e-5 x), compared against a 1e-12 μ floor. The test now passes a partial that matches to 1e-6 relative OR
  whose effect x |∂μ/∂x − fd| is below 1e-9 μ; with it, 13410/13410 pass. A mutant dropping F's partials
  (μ through N_sc,eff) fails 1849 checks, so the test still detects a wrong partial.

| ID | Question | Options | Ruling |
|---|---|---|---|
| KQ8 | KC5's floor for partials below the difference noise | (a) keep the 1e-9 μ effect floor (recommended: the 1e-12 floor was under the FD's own rounding, a test defect rather than a model one); (b) a tighter check by complex-step or higher-precision differences, KC5 open until then | **(a)**, the owner, 2026-11-25 ("choose recommended") |

**Slice 0 record, `klaassen_equilibrium` (2026-11-25).** `ElectronMobility::klaassen_equilibrium` = 3,
`HoleMobility::klaassen_equilibrium` = 4 (KD7); the silicon refusal (KD3); the isothermal node precompute
(`klaassen_at` once per carrier per run) and the electrothermal node table (μ and the total ∂μ/∂T with n_i's
slope; `klaassen_at` per node per evaluation, a cost to revisit if electrothermal runs slow); the deck values,
vocabulary and hints; deck error messages now list all values (KC8). Tests, all Release:
- KC6 met. FD Jacobian, isothermal: 1D heavily doped diode (1e19 / 1e18) worst column 4.5e-9, 2D with field
  mobility 3.0e-9; electrothermal: 1D and 2D-with-field-mobility cases added to the Unit 23 gate, ≤ 1e-6, pass.
  The residual differs from Masetti's at the same state (the model is in use). The Unit 9 diode 0–0.8 V: every
  bias converges, Newton steps per bias equal Masetti's (1, 4, then 5); J is 2.2 % (0.05 V) to 6.2 % (0.5 V)
  below Masetti's. Recorded cause: at 1e17 the Klaassen minority mobilities are below the majority ones
  (electrons 665.6 vs Masetti 721.3, holes 310.7 vs 319.1 cm²/(V s); electron–hole scattering), crossing above
  near 1e18 (313.4 vs 277.0, 184.0 vs 160.3).
- KC8 partly met: the region refusal for either carrier; the deck messages. KC9 partly met: the run record carries
  electrons 3, holes 4 (deck test).
- Binaries run whole: physics, assemble, deck (before the new deck case; the case itself after), app — all pass.
  The full Release and Debug suites have not run (owner's word; no Debug builds by the standing rule).

**KQ8 ruled (2026-11-25):** the owner, "choose recommended": (a), the 1e-9 μ effect floor stands; KC5 met.

| ID | Question | Options | Ruling |
|---|---|---|---|
| KQ9 | `klaassen` with Lombardi surface mobility (KD5 had it chained; found in slice 1 that Lombardi's E_perp slope and the n, p partials would both run through its end mobilities) | (a) refused in this unit, a recorded follow-up, as KD6's others; (b) chained now | **(a)**, ruled as recommended (2026-11-25, standing instruction); KD5's "surface mobility chain" is withdrawn |

**Slice 1 record, `klaassen` (2026-11-25).** `ElectronMobility::klaassen` = 2, `HoleMobility::klaassen` = 3;
`SystemData::klaassen_fluxes` (edge_fluxes with unit factors, then a = (a0 / μ0) μ(n, p, E): harmonic mean, Canali of
it or between materials the mean of each end's Canali); `EdgeFlux::d_other1/2`; the carrier-flux term's pattern and
stamping of the cross entries (4 per edge). Densities below 0 (Newton iterates) continue linearly with the slope at
0, the mobility held at ≥ half its value at 0: the first FD gate run failed (worst 0.5) only because the probe's
1e-7 step took tiny minority densities (9.4e-15 scaled) below 0, where a clamp at 0 halved the difference; the
Jacobian entry itself was right (diagnosed, not a criterion change). Results, Release:
- KC7 met. FD Jacobian with cross entries: 1D heavy 4.5e-9, 2D field mobility 5.7e-8, 1D electrons only 3.1e-9,
  3D 4.3e-8 (≤ 1e-6); the pattern grows by exactly 4 per edge. Unit 9 diode 0–0.8 V: Newton steps equal Masetti's
  at every bias (≤ +2 allowed). High injection (p+ 1e19 / n 1e15, to 0.9 V): base holes reach 1.906e17 cm⁻³
  (191 N_D); `klaassen` against `klaassen_equilibrium` J is +0.4 % at 0.2–0.3 V, −0.26 % at 0.6 V, −11.0 % at 0.8 V
  and −22.8 % at 0.9 V (carrier–carrier scattering); Newton 5–6 both.
- KC8 met: refusals with electrothermal, energy balance, impact ionization, surface mobility (KQ9), quantum
  confinement (each named), `klaassen_equilibrium` allowed with them; non-silicon regions named; deck messages list
  all five / four values.
- KC9 met: run records electrons 2, holes 3 (deck test); a `klaassen` transient (diode pulse, BDF2 adaptive)
  resumes bit for bit at four stop points (Unit 61 harness).
- Run: physics, assemble, app binaries whole; solve `[klaassen]`, `[resume]`; the deck's mobility cases. Full
  suites not run (owner's word); no Debug.

**Slice 2 docs (2026-11-25).** ARCHITECTURE.md: the Unit 66 as-built block at the end of 6.2 and the 11.3 nanoscale
row's mobility mention. Remaining: KC10 (full Release and Debug suites), on the owner's word.

**KC1 met (2026-11-25, owner's word).** NG1 traces of the five traced binaries (assemble, deck, adapt, validation,
solve; Release, MKL) on the branch against `main` 42b4306 (the merge base; its HC1 traces, PARDISO cases from its
MKL re-trace), compared order-independently by case. All five binaries pass. 445 common cases: 443 identical; the
electrothermal FD-gate case keeps every `main` record and adds only its two new Klaassen subcases; the GMRES-ILUT
quasi-static C-V case differs only by attribution: the shared static fixture of `mos_test.cpp` (computed once, by
whichever case runs first) now lands there, not in P3, and the records together are the same. 9 cases are new
(Klaassen 5; PARDISO 4, skipped in the non-MKL `main` run) and P3's key is gone (the attribution above); the 4 PARDISO cases match `main`'s
MKL trace with the same fixture attribution. Run identities (deck traces 65/65 identical) and Jacobian patterns
unchanged.

**KC10 (2026-11-25, owner's word).** Release: build `/W4 /WX` clean, full suite 302/302 (3787 s, ctest -j 4).
Debug not run: the owner stopped it before it started ("dont start the debug process"); Debug is left to CI.

### Unit 67 — Hurkx trap-assisted tunnelling, opt-in (proposed 2026-11-25; on `physics/trap-assisted-tunnelling`)

**Owner request (2026-11-25):** "any physics model?", then "plan for Hurkx trap-assisted tunnelling", "start it in
new worktree" (worktree `nitcad-tat`, branch from `origin/main` 0f49230). Unit 20 recorded trap-assisted tunnelling as
excluded. Without it, SRH lifetimes ignore the field. In a depletion region above about 3e5 V/cm, electrons and holes
tunnel between the bands and the trap (phonon-assisted, through Dirac wells), which shortens the effective lifetimes.
Reverse diode leakage and MOSFET drain/body off-current then come out too low, and their voltage dependence is too
weak.

**Standing conditions.** Opt-in. Unselected, everything is bit-identical to `main`: residuals, Jacobian patterns and
values, run identities, result files, the Unit 40/59 traces and decks. Criteria are fixed before code and never
loosened after a run (TVM2). Builds run `-j 6` and runs use 6 threads. Release suite with `/W4 /WX` runs only on the
owner's word. Debug builds and runs need the owner's word. Commits stay local; no push, PR or merge without the
owner's word. TAQ1–TAQ7 were ruled as recommended under the owner's standing instruction (2026-11-11).

**Sources.**
- (S4) Silvaco ATLAS User's Manual (old edition, `C:\Users\disha\papers\klaassen\atlas.pdf`, not committed),
  "Trap-Assisted Tunneling", eqs 3-71–3-75, Table 3-7 (MASS.TUNNEL 0.25), citing G. A. M. Hurkx, D. B. M. Klaassen,
  M. P. G. Knuvers, IEEE Trans. Electron Devices 39, 331 (1992).
- The 1992 paper itself was NOT read; no copy is in hand. The Charon manual in the same folder implements Schenk's
  band-to-trap model, not Hurkx's, so it is not a source here.
- K below is the WKB exponent of a triangular barrier, (4/3)√(2m) ΔE^{3/2}/(qħF). The PDF extraction of S4 eq 3-74
  is ambiguous. Taken as written below for two reasons:
  - It is the standard exponent.
  - Laplace's method on Γ then recovers Hurkx's widely quoted low-field form exactly: 2√(3π)(F/F_Γ)exp((F/F_Γ)²),
    with F_Γ = √(24 m (kT)³)/(qħ). Derivation: the peak is at u* = (2a/3K)², φ(u*) = 4a³/(27K²) = (F/F_Γ)²,
    and |φ''(u*)| = 9K²/(8a).

**Model.** SRH with a midgap trap, as built (n₁ = p₁ = n_ie). Each lifetime is divided by its field-enhancement
factor (S4 eq 3-71):

  R = (n p − n_ie²) / [ τ_p/(1 + Γ_p) (n + n_ie) + τ_n/(1 + Γ_n) (p + n_ie) ]
  Γ_c(F) = a_c ∫₀¹ exp(a_c u − K_c u^{3/2}) du,   a_c = ΔE_c / kT,   K_c = (4/3) √(2 m_c m₀ ΔE_c³) / (q ħ F)
  ΔE_n = E_c − E_i = kT ln(N_c/n_ie),   ΔE_p = E_i − E_v = kT ln(N_v/n_ie)   (TAQ2)

- F = |E| at the node, the least-squares field of local band-to-band tunnelling (Unit 20).
- F = 0 gives Γ = 0.
- Γ is computed in SI: F is converted from V/cm, ΔE in J, m in kg.
- So a_n and a_p are exactly the node's `log_dos_n_` and `log_dos_p_`. They already carry band-gap narrowing and
  MLDA's band data, so no new tables are needed.
- Silicon at 300 K with m = 0.25: F_Γ ≈ 3.7e5 V/cm.

**Design (TAD1–TAD9).**

| ID | Decision | Reason | Limitation |
|---|---|---|---|
| TAD1 | Physics: new `include/NiTCAD/physics/trap_tunnelling.hpp` and `src/physics/trap_tunnelling.cpp`. They hold `struct TrapTunnellingFactor { double gamma; double d_dF; }` and `[[nodiscard]] TrapTunnellingFactor hurkx_enhancement(double a, double mass_ratio, double T_K, double F_V_per_cm) noexcept`. `d_dF` is in per V/cm. Quadrature is fixed composite Gauss–Legendre: 16 equal panels of 8 points on [0, 1]. dΓ/dF = (K/F)·a ∫ u^{3/2} exp(a u − K u^{3/2}) du on the same nodes. Γ = d_dF = 0 when F = 0 or K > a + 745 (underflow). | Deterministic and bounded. For a peak inside (0, 1), K < 2a/3, so its width is √(8a/9)/K > √(2/a) ≈ 0.18 even at a = 62 (4H-SiC). Each panel then holds a smooth integrand, at most exp(a/16) across it. | Fixed cost of 2 × 128 exponentials per node per evaluation |
| TAD2 | `[models] trap_assisted_tunnelling = off \| hurkx` (default `off`), plus `tat_mass_n` and `tat_mass_p` (default 0.25, S4 Table 3-7). Both masses must be > 0 and finite. One pair for every semiconductor region (TAQ3). | Follows the MLDA mass keys (Unit 47). | No per-material masses |
| TAD3 | Node table: the per-node field rows built for `btbt_local` move into one helper, `build_field_nodes(...)` in drift_diffusion.cpp. It returns the nodes (semiconductor, off the ohmic contacts, with carrier edges) and their `FieldEdge`s and pseudo-inverses. `btbt_local` and TAT share it: `TunnelNode` keeps A and B, and the new `SystemData::TrapNode { node, first, last, inverse[9] }` uses the same `field_edges_`. | One field definition (Unit 20's OLD/NEW). The refactor is under TAC1. | — |
| TAD4 | Term: new `SystemData::Terms::TrapTunnellingSrh` (TermDomain::global; rows n and p; columns n, p of the node and ψ of the node and its edge neighbours). It replaces `Srh` when `trap_assisted_tunnelling = hurkx` and `srh` is on. `Srh`'s predicate becomes `m.srh && !m.tat`. Nodes outside the TrapNode table (contacts, nodes with no carrier edges) get plain SRH, Γ = 0. ψ partials go through `e.ab`/`e.ba` as in LocalTunnelling; the pattern already holds them. | Exact Jacobian; plain SRH stays untouched when off. | — |
| TAD5 | Refused (`invalid_input`, before solving): `hurkx` with `srh = off`, electrothermal, energy balance, density gradient; masses ≤ 0 or non-finite. Allowed with: impact ionization, local and nonlocal band-to-band tunnelling, MLDA, Lombardi, Klaassen, transient, small signal, traces, restart. | Γ depends on T (needs ∂/∂T, a follow-up). Energy balance and density gradient keep their own SRH paths. | Recorded follow-ups |
| TAD6 | Scope is bulk SRH only. Interface SRH (interface_edges.cpp) and interface traps are unchanged (TAQ7). | Hurkx's model is a bulk lifetime factor. | — |
| TAD7 | Outputs: `srh_rate` and RateMechanism::srh report the enhanced rate. No new derived field (TAQ5). | The trap rate is the rate. | Γ not plotted |
| TAD8 | Run records: `models.trap_assisted_tunnelling` = 1 (hurkx), plus `models.tat_mass_n` and `models.tat_mass_p`, all written only when on. Deck vocabulary and help text in `deck_vocabulary.cpp`. ARCHITECTURE.md: an as-built block after Unit 20 and the 11.3 row. | Run identities of existing decks unchanged. | — |
| TAD9 | No new bundled deck. Device checks use the bundled 1D diode deck in tests only (TAC8). | The tree is the owner's. | — |

**Questions (ruled as recommended, 2026-11-25, standing instruction).**

| ID | Question | Options | Ruling |
|---|---|---|---|
| TAQ1 | Form of Γ | (a) Hurkx's integral by fixed composite Gauss–Legendre; (b) the low-field closed form 2√(3π)(F/F_Γ)e^{(F/F_Γ)²} (unbounded at high field); (c) Schenk's model (S1 Charon) | **(a)**; (b) is used only as a check (TAC4) |
| TAQ2 | Tunnelling energy range ΔE | (a) the trap at the intrinsic level, as built: ΔE_n = E_c − E_i, ΔE_p = E_i − E_v from the node's `log_dos`; (b) Hurkx's bias-dependent range (needs the paper); (c) a deck trap level E_t for SRH as a whole | **(a)**; (b) when the paper is in hand, (c) a separate unit |
| TAQ3 | Masses | (a) `[models] tat_mass_n/p`, 0.25 each (S4); (b) per material | **(a)** |
| TAQ4 | Electrothermal, energy balance, density gradient | (a) refused, follow-ups; (b) silent fall-back to plain SRH | **(a)**: no silent model change |
| TAQ5 | Outputs | (a) `srh_rate` enhanced, plus new fields Γ_n and Γ_p; (b) `srh_rate` enhanced only | **(b)**; (a) a follow-up if asked |
| TAQ6 | Published curves | (a) none (paper not in hand); TAC2–TAC8 check limits, the asymptote and device consistency; (b) the owner supplies Hurkx 1992 first | **(a)**, recorded as a limitation |
| TAQ7 | Interface SRH | (a) bulk only; (b) interface recombination enhanced too | **(a)** |

**Criteria (TAC1–TAC11), fixed before code.**

| ID | Criterion |
|---|---|
| TAC1 | Off: every existing test passes unchanged. The Unit 59 NG1 trace set is byte-identical to `main` 0f49230's (built as for KC1: build/ng1 traces). Run identities and Jacobian patterns are unchanged. `btbt_local` cases are included (TAD3 refactor). |
| TAC2 | Γ(0) = 0 and d_dF(0) = 0 exactly. Γ ≥ 0 and strictly increasing on F = 1e3…1e7 V/cm (40 log points) for a ∈ {10, 21.7, 40, 62}, T ∈ {250, 300, 400} K. Finite at 1e8 V/cm, with Γ ≤ a e^a. |
| TAC3 | Quadrature: Γ and d_dF against a 4096-panel × 8-point reference to 1e-12 relative on the TAC2 grid. |
| TAC4 | Laplace check: where u* = (2a/3K)² ≤ 0.25 and (F/F_Γ)² ≥ 5, Γ is within 10% of 2√(3π)(F/F_Γ)exp((F/F_Γ)²). The largest deviation is recorded. Silicon (m 0.25, 300 K): F_Γ = 3.69e5 V/cm to 0.5%. |
| TAC5 | d_dF against central differences of Γ to 1e-6 relative on the TAC2 grid, wherever Γ > 1e-12. |
| TAC6 | Node rate: with Γ_n, Γ_p given, R equals `srh_recombination` with τ_c/(1 + Γ_c) to 1e-14 relative; Γ = 0 reproduces it bitwise. |
| TAC7 | Jacobian: the existing FD-check helper on a reverse-biased (−2 V) diode state, 1D, 2D and 3D, all columns (ψ of the neighbours included), to 1e-6 relative. Also with `btbt_local` and with impact ionization on together. |
| TAC8 | Device: the bundled 1D diode deck, reverse sweep 0 → −5 V in 0.25 V steps, with and without TAT. Requirements: every point converges; Newton iterations per point ≤ the plain-SRH run's + 2; \|I\| with TAT ≥ \|I\| without at every point. The ratio is recorded at −1 and −5 V, and it rises with \|V\|. Forward bias 0.3–0.7 V: the change is recorded. A transient restart continues bit for bit (Unit 61). Small signal at −2 V: the conductance matches the DC dI/dV by central difference to 1e-4 relative. |
| TAC9 | Refusals: each TAD5 refusal names the model or key. The deck error for an unknown value lists `off` and `hurkx`. |
| TAC10 | Run records: the TAD8 keys written when on, absent when off; a two-deck identity check (on vs off) differs. |
| TAC11 | Suites: Release `/W4 /WX` full suite on the owner's word; Debug only on the owner's word. |

**Pre-code corrections (2026-11-25, slice 0, before any code ran; the owner's "go ahead").** I found these
while writing the TAC tests. They replace the plan text above where they differ. They correct defects, and no
criterion was loosened.
- TAD1, quadrature:
  - In u, the low-field integrand peaks within K^{-2/3} of u = 0, next to the u^{3/2} kink, far inside the first
    panel. So the rule is applied after the substitution u = v²: Γ = 2a ∫₀¹ v exp(a v² − K v³) dv, and
    dΓ/dF = (K/F) 2a ∫₀¹ v⁴ exp(…) dv. Both integrands are analytic.
  - The 16 × 8 Gauss–Legendre panels cover [0, v_max] in v. v_max = 1, or the point past the maximum where the
    exponent has fallen 50 below it (bisection).
  - The underflow clause "Γ = 0 when K > a + 745" is wrong and dropped. As F → 0, Γ → a Γ(5/3) K^{−2/3} ∝ F^{2/3}:
    it never underflows, and dΓ/dF is unbounded there. Γ = d_dF = 0 only at F = 0 or where K overflows.
  - So Γ has a cusp (∝ |E|^{2/3}) at zero field. For example, a = 21.7 gives Γ ≈ 0.15 at 1e4 V/cm and ≈ 0.03 at
    1e3 V/cm, so neutral regions with fields of 1e2–1e3 V/cm see a few-percent enhancement. That is the model, and
    TAC8 records the forward-bias change.
- TAC3: the reference is independent of the production rule. It is composite Simpson in v with 2¹⁶ intervals
  (same 1e-12), not a finer copy of the same Gauss nodes.
- TAC4 as written was vacuous. Its region (u* ≤ 0.25 with (F/F_Γ)² ≥ 5) is empty for a < 60, because
  (F/F_Γ)² = a v*²/3 with v* = √u*. Replaced by:
  - (i) the low-field limit: Γ within 1e-3 of a Γ(5/3) K^{−2/3} where a K^{−2/3} = 1e-4, on the TAC2 (a, T) grid;
  - (ii) the Laplace form at v* = 0.7, for a ∈ {40, 62, 100, 200} at 300 K: the relative deviation strictly
    decreasing and ≤ 5% at a = 200, every value recorded;
  - (iii) F_Γ = 3.69e5 V/cm to 0.5% (silicon, m 0.25, 300 K), kept.

**Slice 0 (2026-11-25).** `hurkx_enhancement` with the tests TAC2–TAC5 and `tat extreme field finite`, all
red first against a stub that returns 0. Three cases passed on the stub (zero field, extreme field, and the 1e-300
field), and they guard limits a stub satisfies trivially.
- **The run exposed a quadrature fault, fixed.** With 16 panels, TAC3 missed by up to 1.7e-10 at high field
  (a = 62). There the exponent changes by about 8 across a panel. The panel count is now the larger of 16 and
  ⌈(2a v_max + 3K v_max²) v_max / 4⌉, so φ changes by at most 4 across a panel. The Simpson reference's own error
  there is about 1e-13.
- **TAC5's difference changed.** Central difference with h = 1e-4 F missed by up to 1.6e-6 (a = 62, F ≥ 9e5 V/cm).
  That is the difference's own truncation, about (2(F/F_Γ)² h/F)²/6, because Γ grows like e^{(F/F_Γ)²}. It is now
  the fourth-order 5-point difference, same 1e-6. This is a defect in the test, not the criterion.
- Results:
  - TAC2: met.
  - TAC3: largest relative deviation 5.7e-14.
  - TAC4 (i): met. TAC4 (ii) Laplace deviations: 0.61%, 0.30%, 0.18%, 0.086% at a = 40, 62, 100, 200
    (decreasing, ≤ 5%). TAC4 (iii): met.
  - TAC5: met.
- Physics binary: 139/139 cases.
- Hand mutants, 5/5 caught:
  - d_dF with v³ for v⁴;
  - a 4-point rule;
  - a prefactor dropped;
  - ħ → h;
  - tail 50 → 5.
- The plan's "15 panels" mutant was replaced. With adaptive panels it changes Γ by about 1e-14, so it is
  equivalent under TAC3's 1e-12 and cannot be caught.

**Slice 1 (2026-11-25).** Assembly (`TrapTunnellingSrh`, a global term over all three columns; plain `Srh` when
the model is off), the deck keys, the run record and the vocabularies.
- Results:
  - TAC6 (rate): with each lifetime over 1 + Γ, the term equals SRH to 1e-14 (a zero-current state, unequal
    masses 0.2 and 0.6).
  - TAC7 (Jacobian): the finite-difference gate, worst 3.2e-7; with local band-to-band and impact ionization too.
  - TAC9, TAC10: the deck values, errors, refusals and run record; `models.trap_assisted_tunnelling` and the
    masses are written only when on.
  - Review focus 1 and 2: `tat zero field keeps SRH`, `tat contact nodes plain srh`.
  - TAC1: **met.** NG1 traces of the five traced binaries (assemble, deck, adapt, validation, solve; Release),
    branch against `main` 0f49230, compared order-independently by case. All five pass. 634 common cases:
    632 identical; the two quasi-static C-V cases of `mos_test.cpp` differ only by attribution (the shared
    static fixture lands in whichever runs first), and their pooled records are identical, as in KC1. 9 cases
    are new (the TAT cases).
  - TAC8: **met as amended (below).** The deck tests' diode, reverse 0 → −5 V in 0.25 V steps: every point
    converges, Newton iterations ≤ plain SRH's + 1, |I| with TAT above |I| without at every point, the ratio
    rising at every point: 2.98, 3.50, 4.01, 4.55, 5.11 at −1 … −5 V. Forward: +2.95 % at 0.3 V, +0.42 %,
    +0.060 %, +0.009 %, +0.001 % at 0.7 V. The transient restart continues bit for bit.
- **TAC8 amended by the owner (2026-11-25), on two pre-existing faults the first TAC8 run exposed.**
  1. *The terminal current does not resolve this leakage.* The contact-edge current is a difference of large
     drift and diffusion terms; its own resolution (`terminal_current_resolution`) is about 4e-7 A/cm² on a
     10 nm mesh and 4.5e-5 on a 1 nm mesh, against a generation current of about 5e-9. Read that way the
     reverse current stepped with the mesh and differed 150× between meshes, and its slope was the same
     9.99e-13 S/cm² everywhere. The weighted current (Unit 32) resolves it: 1e-21 resolution, the two meshes
     agreeing to 0.5 %, equal to q∫R. **Ruling:** TAC8 reads the weighted current, and the ratio is checked
     at every point. The small-signal check (no weighted admittance exists) is |Y − G| within Y's own reported
     resolution instead of 1e-4 relative. Recorded as weak: at −2 V the resolution is 1.5e-7 S/cm² against
     G = 8.0e-9, so it constrains little.
  2. *The deck tests' "1D diode" was a compensated resistor.* A profile's window on its depth axis was dropped,
     and the depth axis defaults to x, so in 1D `x=-1:1um` was silently ignored: Na = Nd = 1e17 at every node
     (and the same in `hand_diode`, so deck and hand agreed). **Ruling:** fixed in this unit.
     - `device::evaluate`: a constant profile is bounded by its window on every axis, its depth axis
       included (a constant has no depth shape). `check_profile` refuses a bounded window on the depth axis of
       a Gaussian or erfc (the shape sets it). `deck::profile_text` writes a constant's depth-axis window; the
       app's `profile_extent` reads it. No shipped deck and no app-built profile windows a depth axis
       (2D/3D decks use `axis=y` for x windows), so nothing shipped changes.
     - The deck and app test diodes are now real junctions. Four mechanics tests written against the
       compensated device keep it, bit-identical, as `compensated_device` / `hand_diode(true)`, because their
       tolerances assume its large currents: RG3 (initial checkpoint), the sweep restart, the `[trace]` deck
       (tol_update 1e-12) and RG6 (small-signal restart). The app's 0 V I-V check compares with the point's
       terminal-current resolution instead of exact 0.
     - Doping tests: the 1D step, the refusal, a Gaussian's window unchanged.
- Affected cases run by name (every deck, app and window case in a file using the diode fixtures or an
  x window, hidden cases excluded): all pass. The hidden benchmark "reduced FinFET as layered generate3d
  meshes and is certified" fails P6, and fails the same on `main`'s binary: pre-existing, outside this unit.
  The app completion count is 27 (28 `[models]` keys less srh).
- Hand mutants, 4/4 caught:
  - Γ_n applied to τ_p (the rate test);
  - the other end's ψ columns scaled to 0 (the Jacobian gate; the plan's "dropped" form does not compile
    under `/WX`);
  - the electrothermal refusal removed (the refusals test);
  - the run-record key written when off (the deck run-record test).

**Review focus (not otherwise exercised; each pinned by a test in the slice that owns the code).**
1. A node where F is exactly 0, e.g. a symmetric centre: plain SRH must still be added there. The E_r/F division
   must not run, and LocalTunnelling's `continue` pattern must not be copied. Test `tat zero field keeps SRH`
   (slice 1).
2. Ohmic contact nodes and nodes with no carrier edges: Γ = 0 and plain SRH. Test `tat contact nodes plain srh`
   (slice 1).
3. Extreme field, F ≥ 1e7 V/cm: finite, no NaN, Newton still converges (TAC2 plus test
   `tat extreme field finite`, slice 0).
4. A wide gap (a ≈ 62, 4H-SiC): quadrature accuracy holds (TAC3 grid; slice 0).
5. Small-signal and restart paths read the new ψ columns: the AC Jacobian pattern holds them (TAC8, slice 1).

**Slices and tasks.** Each task follows test → fail → implement → pass → commit (local).

- **Slice 0: physics (TAD1).**
  - Files: create `include/NiTCAD/physics/trap_tunnelling.hpp`, `src/physics/trap_tunnelling.cpp` and
    `tests/physics/trap_tunnelling_test.cpp`; register them in the physics CMake lists.
  - Tests (Catch2, tag `[trap_tunnelling]`): `hurkx zero field`, `hurkx monotone and bounded`,
    `hurkx quadrature reference` (TAC3), `hurkx laplace asymptote` (TAC4), `hurkx derivative` (TAC5) and
    `tat extreme field finite`.
  - Run: `ctest -R trap_tunnelling` in build/release.
  - Hand mutants: drop the u^{3/2} factor in d_dF; use 15 panels; drop the a prefactor; swap ħ for h. Each must
    be caught.
- **Slice 1: assembly, deck, records (TAD2–TAD8).**
  - Files:
    - `include/NiTCAD/assemble/models.hpp`: `bool tat`, `double tat_mass_n`, `double tat_mass_p`.
    - `src/deck/build.cpp` and `src/app/deck_vocabulary.cpp`.
    - `src/assemble/system_data.hpp`: `TrapNode` and `trap_nodes_`.
    - `src/assemble/drift_diffusion.cpp`: `build_field_nodes` and the refusals.
    - `src/assemble/terms.cpp`: `TrapTunnellingSrh` and the predicate.
    - `src/solve/run_record.cpp`.
  - Tests in the assemble and solve test files: TAC6, TAC7, TAC9, TAC10, plus review-focus items 1, 2 and 5.
  - TAC8 runs by case name only (never the whole deck binary).
  - TAC1 trace comparison against main's NG1 traces.
  - Hand mutants: Γ_n applied to τ_p; a dropped ψ column of the neighbour; the refusal for electrothermal
    removed; the run-record key written when off.
- **Slice 2: docs and suites.** ARCHITECTURE.md as-built block and the 11.3 row; this record completed; TAC11 on the
  owner's word.

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
