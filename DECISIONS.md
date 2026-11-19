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

---

## Open units (records kept here in full)

None on `main`. Unit 59 (natural-test norm; approved 2026-11-18, not started) is on `solve/natural-norm`.

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

**Slices (each one local commit):** 26h·1 `axis_box`, drawn, its switch (JD1; JG2, JG6 part); 26h·2 `view_cube`, hit
testing, snapping (JD2, JD3; JG3, JG6 part); 26h·3 named views and smooth moves (JD4, JD5; JG4, JG5, JG8); 26h·4
references, suites, mutants (JG1, JG7, JG9), the ARCHITECTURE entry.

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
