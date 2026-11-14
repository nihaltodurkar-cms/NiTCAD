<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

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

