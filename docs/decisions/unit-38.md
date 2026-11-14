<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

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

