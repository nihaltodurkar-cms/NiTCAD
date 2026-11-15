<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

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
