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
