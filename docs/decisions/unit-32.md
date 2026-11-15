<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

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

