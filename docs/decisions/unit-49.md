<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

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
