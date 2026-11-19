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
