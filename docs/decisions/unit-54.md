<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 54 — 3D structured boundary layers: swept meshes with structured cross-section layers (proposed 2026-10-31; on `mesh/structured-layers-3d`)

**Owner request:** "Implement 3D structured boundary layers (item 5)" (2026-10-31): ARCHITECTURE.md 11.3, FinFET
geometry and meshing, item 5. The branch starts from `mesh/finfet-geometry` `3f08827` (Unit 53, PR #74, not yet
merged), whose complex, extrusion and deck it extends. The standing conditions hold: no change to any existing mesh,
device or result; criteria fixed before code; Release and Debug suites at every commit, `/W4 /WX`; hand mutants;
nothing pushed or merged without authorization.

**Why this unit.** Unit 53's FC9 failed: the generated reduced FinFET (about 2.9e5 isotropic nodes at s = 0.08 nm) is
13.9% low in on-current, 31 mV low in V_T and 28% low in max g_m against the 2× tensor reference (45,375 nodes), while
the coarse tensor mesh (9,367 nodes) is within 1.2%. The 2D MOSFET met the same rule only with a structured layer
(U14). Isotropic tetrahedra at the inversion layer's scale (0.1–0.2 nm normal to the interface) pay that spacing
along the interface and along the fin too; a tensor mesh does not. A FinFET channel needs two anisotropies: (a) fine
normal to the gate stack, coarse along it in the cross-section; (b) fine in the cross-section, coarse along the fin.

**The question that shapes the design (UQ1, below).** Two routes give anisotropic layers in 3D:
- **Route A, swept meshes (recommended).** For an `extrude3d` geometry, mesh the cross-section in 2D (with structured
  layers there) and sweep it along the axis on a 1D line set. The product of a 2D Delaunay triangulation and a 1D grid
  is Delaunay, and its Voronoi dual is the product of the two duals, so the box method is exact by construction
  (argument in UD1). Anisotropy (b) comes free; (a) is U14's 2D layer extended to interfaces and corners. It covers
  both FinFET decks (no void sections) and every extrusion without voids next to filled cells. It does not cover
  general complexes (the API's `PolyhedralGeometry`).
- **Route B, layers inside the 3D Delaunay generator.** A tensor slab on a planar facet or interface piece, embedded in
  Unit 53's conforming refinement, with blocks at corners and junction lines. General, but each slab's circumspheres
  reach past its last row, so the surrounding refinement must be fenced (U14's 2D rule, now with 3D corner and line
  junctions), and termination and P4 at the slab edges are new work. Larger (XL) and riskier; an extruded FinFET
  would still want route A's cross-section layers.

The design below is route A. If the owner picks route B, the design is redone before code.

**Sources.** The Voronoi diagram of a Cartesian product of point sets under the product (Euclidean) metric is the
product of the factors' diagrams (elementary; written out in UD1). U14 and Unit 49 (2D layers), FD2 (the extrusion),
FD7–FD10 (interfaces, patches and their checks), M7 (mirrored pairs), M3 (the box method's conditions).

**Design:**

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| UD1 | **The swept mesh** (`mesh::generate_swept(ExtrudedGeometry, SweptOptions)`, API and deck): the cross-section mesh M₂ from the 2D generator; axis lines t₀ < … < t_m; nodes (p, t_j) for every 2D node p and line t_j inside the domain; each prism (triangle × [t_j, t_{j+1}]) split into three tetrahedra by the global-index rule (every quad face cut on the diagonal through its smallest node index), so neighbouring prisms conform. The tetrahedra go through `make_dual_mesh` and its checks (P0–P3, P5) like every generated mesh. **Argument:** a prism's six vertices lie on the sphere centred at (c, t_{j+½}) with radius √(R² + Δ²/4), c and R the triangle's circumcircle; a node (q, t_k) of any other prism is outside or on it, since \|q − c\| ≥ R (M₂ Delaunay) and \|t_k − t_{j+½}\| ≥ Δ/2. So the tetrahedralization is Delaunay with ties exactly inside each prism, and its dual is V₂(p) × [t_{j−½}, t_{j+½}]: volumes A₂(p)·Δt_j, in-plane couplings (ℓ₂/d₂)·Δt_j, axial couplings A₂(p)/Δ. | Proven (above); P1 holds with the ties decided exactly by the perturbed insphere. Checked: UC3 compares every volume and coupling with the product formula. | Fine in the cross-section, coarse along the fin, the box method exact; no new 3D refinement. | Extrusions only. Tie diagonals give zero-area dual faces (rounding judged by P3's τ, as with 2D layers). |
| UD2 | **The cross-section** is meshed with one 2D polygon per cell. Cells with the same region in every section (the same column of regions) share a 2D region, so their common edges are internal; every other cell edge is a 2D interface (M7 pairs at ±s/2; the 2D rules: 6 s, 4 s, 60°). Every lateral interface of any section is then straddled in M₂. | M7, U12, U15 (2D; proven away from features). | Every section's lateral interfaces are among the straddled cell edges. | Some nodes straddle a line that is an interface only in another section. |
| UD3 | **Axis lines:** given by the user (lengths); the two end stations are lines exactly; every interior station where some cell changes region is **straddled**: lines at station ± w/2 (w given per station, as the tensor deck's `straddle`), no line strictly between them. A station where nothing changes needs nothing. A node's 3D region is its 2D cell's region in the section holding its line. | M7 in 1D: the dual face between the two lines is the station plane. | The gate edges (spacer \| HfO₂ at x = 30 and 55 nm) are exact interfaces. | — |
| UD4 | **Refused in this unit** (`invalid_input`, naming the cell and section): a void (cell, section) next to a filled one across a station; that boundary would need nodes on the station for some cells and straddles for others. A cell void along the whole axis is simply not meshed. | — | Both FinFET decks have none. | Raised source/drain ending mid-axis next to void needs Unit 53's generator or a later unit. |
| UD5 | **2D structured layers on interfaces** (U14 extended; `BoundaryLayer` gains signed lines): a layer is a base segment (a boundary segment, as now, or an interface segment: a chain of collinear cell edges), positions `along` it, and **signed normal lines** (offsets from the base, positive to its left, negative to its right). Every interface parallel to the base inside the layer's span is straddled by two lines at ± s_i/2 about it with no line between (s_i a layer option, default s); a boundary parallel to the base inside the span is a line exactly; the lattice lies inside the domain. Cell edges crossing the layer (normal to the base) meet it only between two straddling columns: `along` straddles each crossing likewise. No M7 pairs are placed on an interface piece a layer covers; the lattice replaces them there. | As U14: the lattice's rectangles are Delaunay ties, its dual the tensor box method; a straddled interface lies on the dual face midway (P4a exact). | The gate stack (0.5 nm SiO₂, 1.5 nm HfO₂) and the inversion layer under it in one lattice, as the tensor deck's `straddle` lines do. | The user writes the lines (as for the tensor deck); `linspace` and straddles help. |
| UD6 | **Corners where two layers meet** (the fin's top corners, its foot): at a **right angle** the corner square is the tensor product of the two layers' normal lines (a corner block), and each layer's `along` starts with the other's lines, so the three lattices are one tensor lattice: exact. At **other angles in [60°, 120°]** (the tapered fin's 95° top, 85° foot) a parallelogram block: nodes where lines parallel to one layer cross lines parallel to the other, Delaunay-triangulated; its interface bends are straddled to s by construction, but its dual faces tilt by the angle's departure from 90°, so P4b judges it as a corner zone. Below 60° or above 120°: refused. | Right angle: proven (tensor). Other angles: heuristic, checked by P4b (corner deviation ≤ s, measured and pinned). | The channel wraps the fin's corners; the 2D MOSFET had none. | Tapered corners accurate to s, like Unit 53's feature zones. |
| UD7 | **Between layers and the rest:** refinement never inserts into a layer's or block's region and does not judge triangles of lattice nodes alone (U14); the last row should be graded to the neighbouring h (reported, not enforced: a coarse last row makes slivers, as Unit 49 found). Boundary fill between lattice columns keeps Unit 49's rule. | U14, Unit 49. | Reuses the tested 2D machinery. | — |
| UD8 | **Interfaces and patches in 3D:** the product mesh is certified with Unit 53's checks (P4 over the 3D pieces, P6 over patch edges). Lateral pieces are exact wherever M₂ is (UD2, UD5), station pieces exact (UD3); where a lateral piece meets a station piece (a junction line) each dual face lies within one section, so the product of the 2D and 1D straddles is exact there too, measured. A patch over part of the axis (the gate, 30–55 nm on HfO₂'s outer faces) ends at a straddled station: FD8 (a), exact. A patch ending at a station that is not straddled is refused. | Proven by the product structure; checked by P4/P6 as for every multi-region mesh. | The certificate `Device::create` requires, unchanged. | — |
| UD9 | **Options** (`SweptOptions`): the 2D options (h, per-region h, boxes, s, minimum angle, `max_nodes` counted over the 3D nodes) plus the axis `lines`, per-station `straddle` widths, and the layers (UD5), corner blocks inferred where two layers' bases meet. | — | — | — |
| UD10 | **Deck:** in `[mesh] type = extrude3d`, `sweep = on` selects UD1 (off: Unit 53's generator, unchanged); `lines = <lengths>` (`linspace(...)` and `+` as in the tensor deck); `straddle = <station> <width>`; `layer = <y z> <y z> along=<lengths> lines=<signed lengths> [straddle=<length>]`, the endpoints exact cross-section points on cell edges. `nitcad check` reports nodes, prisms, tetrahedra and layers. | As FD13: exact inputs, nothing computed by intersection. | — | — |
| UD11 | **io, run record, app:** a swept mesh is a generated 3D mesh (simplices, regions, patches, certificate) for files, run records (the generation digest also covers lines, straddles and layers) and the application; FD14 applies unchanged. | — | — | — |
| UD12 | **Existing paths:** the 2D generator without layers on interfaces, Unit 53's generator, the convex generator and the tensor path stay bit-identical (UC1). | — | — | — |
| UD13 | **Cost:** M₂'s cost plus a linear product; budget 1e5 nodes ≤ 10 s Release (M₂ about 1e3–1e4 nodes). | — | — | — |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| UC1 | Nothing changes: every existing test and pin passes; every mesh of the 2D generator (boundary layers included), Unit 53's generator, the convex generator and the tensor path keeps its graph, region and certificate digests; every existing deck's run identity and device file byte-identical (`sweep` off). |
| UC2 | Refusals: every UD3, UD4, UD5, UD6, UD8 and UD10 rule refused with its code, index or deck line, and message (one test each). |
| UC3 | Product exactness: on swept meshes of a rectangle, a tapered-fin cross-section and an L-shape, every node volume equals A₂(p)·Δt_j and every coupling the product formula within 1e-12 relative; P0–P3, P5; two runs bit-identical. |
| UC4 | Linear field: u = a·x_k for each axis, Dirichlet on the faces normal to x_k, Neumann elsewhere: conjugate gradients on the swept dual give u at every node within 1e-12 of its range. |
| UC5 | Uniform-along-the-axis equality: a gated tapered fin, uniform along its axis, swept from a 2D mesh: the 3D equilibrium potential on every line equals the 2D solve on M₂ within 1e-9 relative, and the gate charge per unit length within 1e-9. |
| UC6 | 2D layers on interfaces (UD5, UD6): a two-material strip with a layer straddling the interface; a three-layer stack (Si, 0.5 nm, 1.5 nm, boundary); a rectangular fin with right-angle corner blocks; a tapered fin with 95°/85° blocks. Every lattice node bit-identical to its formula; the rectangular parts' dual equal to the tensor box method within 1e-12; P4 certified: flat deviations within P4a, corner deviations ≤ s, measured and pinned. |
| UC7 | 2D accuracy: FC8's cross-section (the flat gate) with layers: discrete H1 potential order ≥ 0.9, L2 ≥ 1.5, gate charge order ≥ 0.9 (Unit 16's thresholds) against a 4× finer layered reference; reported beside Unit 53's FC8 measurement. |
| UC8 | **The reduced FinFET (FC9 re-run; the unit's gate):** the swept deck with layers has at most the coarse tensor mesh's 9,367 nodes; against the 2× tensor reference (45,375 nodes, already measured), I_D at every V_G above the current resolution, V_T, max g_m and the swing satisfy \|Q_gen − Q_ref\| ≤ 1.5 \|Q_tensor − Q_ref\| + f (FC9's rule, unchanged); every point converges. Risk stated before the run: if it fails, it stays `[!shouldfail]` with its measurement for an owner ruling; the threshold is never loosened. |
| UC9 | The shaped FinFET (tapered, wrapped stack, spacers) as a swept deck with layers: certified; the 0–0.8 V sweep converges at every point; nodes, mesh and solve time reported beside Unit 53's FC10 run (288,729 nodes, 3,099 s), and its on-current, V_T and swing compared with that run. No accuracy gate (no reference). |
| UC10 | Deck, io, run record, app: the vocabulary checked by `nitcad check` with line numbers; a swept mesh round-trips with patches, regions and certificate; the run identity covers lines, straddles and layers; `nitcad plot` and the application show it. |
| UC11 | Cost: mesh time at about 1e4 and 1e5 nodes against UD13. |
| UC12 | Hand mutants: at least one per design row (UD1–UD13), each survivor closed by a test or argued equivalent. |
| UC13 | Release and Debug suites pass at every commit, `/W4 /WX`; the Debug suite grows by ≤ 60 s (Unit 53 overran its budget); larger cases go to `[.benchmark]`. |

**Slices** (each a local commit with both suites green):
1. The swept mesh without layers: product, tetrahedra, regions, straddled stations, certificate (UD1–UD4, UD8, UD12;
   UC1, UC2 part, UC3–UC5).
2. 2D layers on interfaces with signed lines and straddles (UD5, UD7; UC6 part).
3. Corner blocks, right-angle and oblique (UD6; UC6, UC7).
4. Deck, io, run record, application; the two FinFET decks as swept decks (UD9–UD11; UC10).
5. Gates and cost: UC8, UC9, UC11, mutants UC12, ARCHITECTURE 6.11 and 11.3.

**Size:** L to XL (5 slices; slices 2 and 3 are the bulk, in `generate2d.cpp`).

**Owner questions (UQ):**
- UQ1: route A (swept meshes for extrusions, structured layers in the cross-section) as designed, or route B (layers
  inside the general 3D generator)?
- UQ2: approve UD1–UD13 and UC1–UC13 as written?
- UQ3: UD6's oblique corner blocks (heuristic, judged by P4b as corner zones) in this unit, or right-angle corners only
  (a tapered fin's layers then stop short of its corners, which get ordinary interface pairs)?
- UQ4: UD4's refusal of voids next to filled cells across a station: accepted for this unit?

**Owner answers (2026-10-31), design and criteria approved:** UQ1 route A (swept meshes); UQ2 UD1–UD13 and UC1–UC13
approved as written, start slice 1; UQ3 oblique corner blocks in this unit; UQ4 voids next to filled cells refused.

**Slice 1 (2026-10-31): the swept mesh without layers** (`mesh::generate_swept`, `generation_digest(ExtrudedGeometry,
SweptOptions)` in `generate.hpp`; `src/mesh/generate_swept.cpp`; `tests/mesh/generate_swept_test.cpp`, 4 cases;
UC5 in `tests/solve/complex_convergence_test.cpp`):
- UD1 as designed: M₂ from `generate_mesh(Geometry2D)`, node (p, j) = j·n₂ + p, three tetrahedra per prism by the
  global-index rule, lateral boundary quads cut on the same diagonal, then `make_dual_mesh` and its checks. The
  extrusion is validated by `extrude` and indexed by `ComplexIndex`; boundary facets' patches are found again by
  their exact vertex sets, so the patch rule stays `extrude`'s alone.
- UD2: a 2D region per distinct (regions per section, both end patches, the set of lateral patches per section).
- UD3/UD4 as designed. **Added, for review:** a station where a lateral patch changes may also be an axis line
  (the patch edge is then a mesh line, exact, no straddle); a lateral patch change at a cross-section vertex where
  the boundary is straight and one 2D region continues is refused (M₂ has no node there).
- UD8: the certificate covers the complex's interface facets as pieces at resolution s; P6 with every edge between
  boundary facets of different patches declared and conservative zones at their ends. **Added, for review:** a swept
  mesh is also certified when it has per-corner patches but no interface (two 2D regions of one 3D region meeting
  on a patched face), per FD10/SQ3.
- Measured: P4 on the gated fin (two straddled stations, three regions, a partial gate patch): flat 7.3e-4 of the
  tolerance, corner deviation 2.2e-16 (the product is exact at junction lines too, as UD8 argues); P6 flat 1.3e-3;
  the gate patch's area exact to 3e-15. UC3 (rectangle, L-shape, tapered fin with an oxide wrap): every volume and
  coupling equal to the product within 1e-12, prism diagonals zero. UC4: every axis within 1e-12. **UC5:** the
  flat-gate fin of FC8, 796 2D nodes × 3 lines: potential equal to the 2D solve within 1.1e-16 V of 0.71 V, gate
  charge per length within 3.9e-16.
- Suites: Release 128/128, Debug 117/117. Debug growth: `mesh.generate_swept` 31.4 s, UC5 about 8.5 s (UC13 budget 60 s).

**Slice 2 (2026-10-31): 2D layers on interfaces with signed lines and straddles** (`BoundaryLayer::lines` and
`straddle` in `generate.hpp`; `read_signed_layer` and `finish_interfaces` in `src/mesh/generate2d.cpp`;
`tests/mesh/generate2d_layers_test.cpp`, 4 cases):
- UD5 as designed. A layer gives `lines` (signed offsets, positive to the left of a → b) in place of `offsets`. The
  generator adds the straddles itself, as UD3 does on the axis: lines at d ± s_i/2 about every interface parallel to
  the base inside the lattice, and columns at t ± s_i/2 about every interface across it (strictly inside the first and
  last columns). A given line or column inside a straddle is refused, and so is a parallel interface without a line
  beyond its straddle on each side. Interfaces meeting the lattice at any other angle are refused until slice 3. A
  covered chain gets no pairs there. Boundaries parallel to the base inside the lattice must be lines, and boundaries
  across it columns. Their nodes are pinned to the edge: a polygon vertex or patch end exactly; on an axis-aligned
  edge only the constant coordinate is set; on an oblique edge the node is projected.
- UD7 as designed: refinement stays out of the rectangle, triangles of layer nodes only are exempt, and there is no
  boundary fill between lattice nodes.
- U14's node order is kept: the pairs are placed after the layers are validated (`finish_interfaces`), and both kinds
  of layer add their nodes after the pairs, as U14 did. The `offsets` path is unchanged, and every existing mesh test
  and pin passes (UC1).
- **Added, for review:**
  (a) The base is any segment covered by boundary and interface edges on its line. A crossing interface splits the
  boundary edge at its junction, so "inside one boundary edge" would refuse a boundary layer that an interface
  crosses.
  (b) In signed mode n is always the left normal of a → b, so a boundary base must have the domain on its left, or
  use negative lines.
  (c) The uncovered part of a covered chain gets pairs graded 1.5 s from the layer's edge and 1.5 s, 3 s from a
  feature, and is refused if it is too short for them. The 6 s rule still applies to whole chains.
  (d) The 4 s separation rules are measured between uncovered parts only. They now run after the layers, so for an
  input that breaks both a separation rule and a layer rule, the layer error comes first.
  (e) A junction inside a signed layer gets no boundary straddle: the lattice's columns or lines straddle it.
  (f) A patch end inside the lattice must lie on a line (it becomes a column, outside every straddle). A domain vertex
  on a lattice line must lie on a column (or be a straddled junction), and one on a lattice column on a line.
  (g) The "too short" refusal is a chain's rule and carries no layer index.
  **Owner ruling (2026-10-31): "choose recommended": (a)-(g) accepted as written.**
- **Measured (UC6 part):** a strip with its interface straddled by one lattice covering the square: the mesh is the
  lattice alone (40 nodes, no pairs), every node and coupling equals the tensor box method within 1e-12, the
  rectangles' diagonals have zero coupling, P4 flat deviation 0. The three-layer stack (Si, 0.5, 1.5 to the
  boundary, s = 0.1): both interfaces straddled, the top boundary a line; P4 flat 0, corner 0. A boundary layer
  crossed by an interface: columns at 0.5 ± s/2, no junction straddle, the pairs resuming 1.5 s above the lattice;
  P4 flat 0, corner 0.
- **Exactness holds where the lattice is Delaunay whatever lies outside** (UD7's grading, reported and not enforced).
  A row matches the tensor values to 1e-12 when every rectangle next to it has its circumdisc inside the lattice, or
  the lattice ends on the boundary on that side. In the stack (columns 2.5, lowest rows 1.0 and 0.8 apart, h = 2),
  the lowest row's rectangles reach the refined Si: rows 0–1 differ by up to 9% in volume, and rows 2 and up are
  exact. The test checks exactly the rows the circumdisc rule admits and pins how many there are.
- Digest: `generation_digest(ExtrudedGeometry, SweptOptions)` also covers `lines` and `straddle`. Slice 1's swept
  digests change; none was pinned.
- Refusals (UC2 part, one test): both offsets and lines, a negative straddle, a base off the edges, fewer than two
  lines, unordered lines or columns, `along` not reaching the end, a line inside a straddle, no line beyond a
  straddle, the protected band, the lattice leaving the domain, a boundary inside the lattice that is not a line, a
  column inside a crossing's straddle, a crossing within its straddle of the first column, an uncovered part too
  short, a domain vertex between columns, an oblique interface, a patch end between lines.
- Suites: Release 129/129, Debug 118/118. Debug growth: `mesh.generate2d_layers` 0.22 s (UC13 budget 60 s; 40 s used by slice 1).

**Slice 3 (2026-10-31): corner blocks, right-angle and oblique** (`src/mesh/generate2d.cpp`: the signed layers are
prepared in four passes, `parse_signed`, `link_corners`, `straddle_signed`, `finish_signed`; `Layer` gains a region
(`tmin`, `tmax`, `cuts`, `polygon`); the header documents the corner rules; 3 cases in
`tests/mesh/generate2d_layers_test.cpp`; UC7 in `tests/solve/complex_convergence_test.cpp`):
- **Changed from UD6 as written, for review (the oblique block).** UD6 put a parallelogram of crossing lines at an
  oblique corner. Its pairs across each interface are offset along it by s tan(angle - 90 degrees) as far from the
  corner as the other layer's lines reach (2-3 nm), far outside P4's corner zone (2 s from a piece's end), so P4a
  would refuse every such mesh. Built instead, one rule for every angle in [60, 120]: a **mitred corner**. Each
  layer's region is cut at the corner's bisector and reaches past its base end up to it; there its columns are the
  other layer's lines (with their straddles) where they cross the bisector, t = t_c + kappa o, kappa =
  +-(1 + cos)/sin of the angle; each layer keeps its nodes on its own side (the lower-indexed layer those on the
  bisector). Its columns stay perpendicular to its base, so its pairs are mirrored right up to the bisector, and only
  the bends near it are judged by P4b. At a right angle kappa is +-1 (exactly, for axis-aligned bases) and the two
  lattices and the corner square are UD6's single tensor lattice, as UD6 prescribed.
- **Added, for review:**
  (a) The two layers at a corner must span the same offsets (equal first and last lines, signs matched across the
  bisector), so the cut regions close on it with no gap or overlap.
  (b) A layer's given columns must stop before the corner (refused strictly between its innermost corner column and
  the base end); the end's own column, which `along` must still contain, is dropped.
  (c) Refused: three layers at one end, the two corners of one layer overlapping, an interface across the layer or a
  patch end inside a corner.
  (d) Triangles of signed layers' nodes alone are exempt from refinement only where their centroid lies in a lattice
  (U14 exempts all of them). A gap enclosed by lattices, as the fin's core between its three layers, is otherwise
  never refined. U14's layers (`offsets`) keep U14's rule (UC1).
  (e) A signed layer's base is on its edges' line within 1e-12 of its length (it was exact): a base end inside an
  oblique edge is a rounded point.
  (f) Layers must not overlap or touch, tested on their regions (separating axes; it was a vertex test that missed two
  crossing lattices); a corner's two layers touch along its bisector.
  (g) The signed layers' nodes now follow every U14 layer's (they are prepared after them).
- **Measured (UC6):** the rectangular fin (top and both sides, two right-angle corners, three interfaces, different
  lines on top and sides): every node of the frame at its formula bit for bit (a corner node at its owning layer's),
  no other node inside the frame, 348 top-band and 179 side-band nodes equal to the tensor box method within 1e-12
  (rectangle diagonals zero), P4 flat 2.0e-3 of the tolerance, corner 1.8e-15. The tapered fin (5 degrees, substrate,
  floor and wrap of uniform thickness, five layers, four 95-degree corners at the top and the foot, no pair left):
  flat 6.6e-3, corner 2.7e-15 with one set of lines, and 2.7e-15 too with other lines on the sides (the straddles about
  each bend are the same on both layers).
- **UC7 fails its H1 threshold (for an owner ruling).** FC8's flat-gate cross-section, layers on the fin's top and
  sides (from y = 14 nm; 93.6-degree corners; 2 nm into the silicon to 4 nm out), lines and columns divided by k,
  against the layered mesh at k = 16 (319,605 nodes, 199 s): H1 orders 1.49, **0.86**; L2 2.28, 1.63; gate charge
  1.83, 1.59 (Unit 53's unlayered 3D FC8: H1 0.43, -0.04). The corner deviation is 3.4e-4 nm at every level (the wrap
  is 0.998 nm thick on the sides, 1 nm on top, so the bends sit just off the bisector; well within s). By zone, the H1
  error squared converges at order 1.10 in the corner blocks, 0.99 in the rest of the lattices and 1.01 at the
  (unlayered) fin feet; the miss is the HfO2's outer bends at (3.375, 10) and (16.625, 10) nm, re-entrant corners
  with a permittivity contrast (22 against 3.9) that no layer covers, at order 0.66 (and the region beyond them).
  The case is `[!shouldfail]`, its thresholds unchanged.
- Refusals (UC2 part, one test): the spans differing, a given column inside a corner at either end, three layers at
  one end, a layer's two corners overlapping, a 135-degree corner, an interface across a layer inside an oblique
  corner, a patch end inside a corner.
- **Owner ruling (2026-10-31): "choose recommended":** the mitred corner accepted in place of UD6's parallelogram;
  additions (a)-(g) accepted as written; UC7's H1 miss accepted as a limitation (the unlayered HfO2 bends), the case
  staying `[!shouldfail]` with its thresholds unchanged.
- Suites: Release 129/129, Debug 118/118 (the Debug run found a test polygon missing a shared vertex in one new
  refusal case; fixed, and the mesh suites re-run 17/17 in both builds). UC7 re-run: "failed as expected", the same
  numbers. Debug growth: `mesh.generate2d_layers` 0.22 s to 2.5 s (UC13 budget 60 s; about 42 s used).

**Slice 4 (2026-10-31): deck, io, run record, application; the two FinFET decks as swept decks** (UD9-UD11; UC10, UC2
part). `[mesh] type = extrude3d` takes `sweep = on | off` (default off, Unit 53's generator unchanged); with it,
`lines` (an axis, as a tensor deck's, required), `straddle = <station> <width>` (repeatable), `layer = <u v>, <u v>
along=<axis> lines=<axis> [straddle=<length>]` (repeatable), `min_angle`, `max_size`, `interface_resolution`,
`max_nodes` and `box` (`src/deck/build.cpp`: `option_swept`, `layer`; the axis parser now takes a line and a text, so
`lines`, `along=` and `lines=` share it; `include/NiTCAD/deck/build.hpp` documents the keys). `nitcad check` prints
the swept counts (axis lines, prisms, tetrahedra, layers, corner blocks, the cross-section's smallest angle) and the
generation digest; the run record, the mesh file and a file mesh read back name the generator as for Unit 53 (FD14
unchanged). The application's vocabulary knows the four keys (and `min_angle` for extrude3d); its mesh menu offers
`sweep` for an extrude3d deck, and `min_angle` instead of `radius_edge` once it is on. `nitcad plot` draws the swept
meshes (mesh, slice and quality views checked by eye on both FinFET decks).
- **Added, for review:**
  (a) **A generator defect fixed (slice 2's code):** a layer whose base line holds two or more collinear interface
  chains (three regions meeting along one line: the fin, the interfacial layer and the HfO2 on the buried oxide's top)
  placed one straddle pair per chain and refused the duplicates as "two interface straddles overlap or touch". One
  pair at the shared offset straddles every chain on the line exactly, so equal straddle centres (within the
  layer's tolerance) are now placed once, for lines and for columns. Only inputs refused before change (UC1). Test:
  three regions meeting at a T on the layer's line give the tensor lattice exactly (flat 0, corner within 1e-12);
  with the merge disabled the test fails with the old refusal (mutant killed).
  (b) P4a's failure message gives the edge's midpoint, as P4b's does (it found the cause in the shaped deck, below).
  (c) `MeshReport` gains `prisms`, `axis_lines`, `layers`, `corner_blocks` (swept meshes only; corner blocks counted
  as end points two signed layers share). The mesh file's report carries them as an optional key (5), written only
  when `prisms` is nonzero, so every other generated file is byte-identical.
  (d) Under `sweep = on`, `box` is the cross-section's (u0 v0 u1 v1 h, as generate2d's) and `max_nodes` counts the 3D
  nodes (the 2D cap keeps its default); `radius_edge` is refused there, and `lines`, `straddle`, `layer` and
  `min_angle` are refused without it, each with its line.
  (e) The run record's mesh type stays `extrude3d` for both generators; the digest tells them apart (the swept
  digest is tagged and covers the lines, straddles and layers).
  (f) `layer`'s value is split at its options' names (`along=`, `lines=`, `straddle=`, at parenthesis depth 0): the
  base's two points before the first, each value running to the next name. Unknown, repeated, empty or missing
  options are refused with the line.
  (g) `along` must still end at the base's length within 1e-12 (slice 2's rule); a tapered side is sqrt(904) nm long,
  which the shaped deck writes as `30.0665927567458nm`. A symbolic end (say `end`) would be kinder; not added.
- **The decks.** `decks/finfet_3d_swept.deck` (the reduced fin, UC8's mesh) and `decks/finfet_3d_shaped_swept.deck`
  (the tapered fin, UC9's), the devices, regions, patches, doping and analysis of Unit 53's two decks; lines 2.5 nm
  apart along the fin, the gate edges straddled by 0.5 nm; h = 4 nm in the cross-section; s = 0.08 nm (the 6 s rule
  applies to every interface piece, covered or not, and the interfacial layer's piece on the box is 0.5 nm).
  - Reduced: four layers wrap the fin (both sides, the top, the foot on the buried oxide), right-angle corners, lines
    from 3 nm into the silicon to the gate (2 nm out; 2 nm into the box at the foot): every interface of the
    cross-section lies inside the frame. 50,061 nodes (1,353 per cross-section x 37 lines), 3.5 s in Release; P4
    flat deviation 0, corner 1.1e-14 cm, P6 present. **Five times UC8's budget (9,367 nodes)**: the four corner
    blocks alone hold 4 x 13 x 13 nodes per section. Slice 5 measures UC8 and tunes the lines first; the risk stated
    in UC8 stands.
  - Shaped: the interfacial layer's outer face is not parallel to the fin's side (0.5 nm wide horizontally at the
    foot, 0.533 nm at the top), so no layer line can follow it: the layers run from 2 nm into the silicon to 0.1 nm
    out and the interfacial | HfO2 interface keeps the 2D generator's mirrored pairs; oblique corners (93.8 degrees at
    the top, 86.2 at the foot). 71,188 nodes, 5.3 s, certified with P6 (Unit 53's generator: 288,729 nodes, about
    20 minutes).
  - **Found while building them (UD7):** with columns 3 nm apart the shaped deck failed P4a: the lattice's thin
    rectangles (rows 0.1-0.6 nm apart) have circumcircles reaching 1.5 nm, past its last row to the pairs 0.36 nm
    beyond, so Delaunay dropped the lattice's edges there and an edge ran from a pair straight across it. UD7's
    grading has to hold along the layer too: columns no farther apart than the gap to the nearest node outside.
    Both decks use 0.5-1 nm columns. UD7 reports rather than enforces this, as written; nothing checks it yet.
- Tests: the deck builds the generator's mesh from the keys exactly (the same digest, graph and regions as
  `generate_swept` on hand-built inputs; each layer line a row of nodes; a cross-section box refines; sweep = off is
  Unit 53's), every mistake with its line (22 cases), the run record, mesh file and file mesh round trip, the run
  identity changed by lines, a layer's straddle and its lines, not by max_nodes; io: a swept mesh with a partial gate
  patch round-trips (corner patches, P6, swept counts, digest), the swept digest changes with each of 24 inputs and
  is not the complex's; app: the vocabulary and menu (WG7 with a swept resistor); both FinFET decks build and
  certify (8.6 s in Release; hidden, below).
- Suites: Release 129/129 (845 s), Debug 118/118 (1,086 s). The FinFET decks' test was visible in that Debug run and
  took `deck.build` from about 10 s to 285 s, so it is `[.benchmark]` (run by name); re-run, `deck.build` takes 15.5 s
  in Debug and 0.7 s in Release. Debug growth this slice: about 5.5 s in `deck.build` and 7 s in `app.deck_sections`
  (WG7's swept resistor); Unit 54 at about 55 s of UC13's 60 s. UC7 re-run in Release after the straddle merge:
  "failed as expected", the same numbers (H1 1.49, 0.86; L2 2.28, 1.63; gate charge 1.83, 1.59).

**Slice 5 (2026-10-31): gates and cost** (UC8, UC9, UC11, UC12; `tests/deck/complex_finfet_test.cpp` `[.uc8]`,
`[.uc9]`, FC9's comparison factored into `compare_to_tensor`, unchanged; `tests/mesh/generate_swept_test.cpp`
`[.benchmark]` UC11):
- **The reduced swept deck tuned to UC8's budget before any solve.** `decks/finfet_3d_swept.deck` takes the coarse
  tensor mesh's lines where they cross the frame (1 nm normal lines; the stack's interfaces straddled by 0.2 nm, the
  foot's by 0.5 nm; columns about 5 nm apart on the sides and top), the rest of the cross-section unstructured at up
  to 8 nm, and the coarse tensor mesh's axis except 10 nm (not 5 nm) spacing in the 1e20 source and drain beyond 20
  and 65 nm: 332 nodes per cross-section × 25 lines = **8,300 nodes** (budget 9,367), certified, generation digest
  `3b758bc4fbfd407a`. Slice 4's lines (50,061 nodes) are gone.
- **UC8 fails (the risk stated before the run); accepted as a limitation by the owner (2026-10-31: "UC8 accept the limitation").** Release, PARDISO on 6 threads; every run
  converged (17/17). Coarse tensor 9,367 nodes (sweep 15 s), reference 45,375 (491 s), swept 8,300 (build 1 s,
  sweep 17 s). **Within their bounds:** V_T 0.2275 against 0.2310 V (bound 6.3 mV, error 3.4 mV), max g_m 7.42e-5
  against 7.63e-5 A/V (−2.7%; the tensor −2.1%), the swing 64.63 against 64.89 mV/dec; I_D at 12 of 17 points,
  0.8 V among them (1.916e-5 against 1.897e-5 A, +1.0%; the tensor +1.2%). **Outside:** I_D at 0, 0.05 and 0.1 V
  (−3.4%, −2.7%, −1.7%; the tensor −1.8%, −1.1%, −0.14%) and at 0.4 and 0.5 V (−1.1%, −0.6%; the tensor −0.6%,
  −0.3%). At 0.1, 0.4 and 0.5 V the tensor's own error is near a sign change, so the bound there is small. Against
  FC9's generated mesh (297,758 nodes: −13.9% in on-current, V_T 31 mV low, max g_m −28%), the swept mesh has
  1/36 of the nodes and is about as accurate as the coarse tensor mesh, but not uniformly closer, as UC8 requires.
- **One retry, stated before it was run:** the fin's core (between the lattices) unstructured at 8 nm, where the
  tensor mesh has lines about 2 nm apart, while subthreshold current flows through the whole fin. A cross-section
  box at 2 nm over the core (9,300 nodes) made it worse: I_D outside its bound at 9 of 17 points, subthreshold
  now +6.7% (0 V) high, on-current +0.2%; V_T, max g_m and the swing still within. The subthreshold current moves
  by several percent either way with the core's resolution, so neither mesh resolves it at this size. No further
  tuning against the reference; the deck is the first one, the case `[!shouldfail]` with its rule unchanged.
  **Owner ruling (2026-10-31): "UC8 accept the limitation".** The swept reduced deck stays as measured (8,300
  nodes), the case stays `[.uc8][!shouldfail]` with FC9's rule unchanged; the subthreshold and near-sign-change I_D
  misses are a recorded limitation of this mesh size.
- **UC9:** **met.** `decks/finfet_3d_shaped_swept.deck` (`[.uc9]`, Release, PARDISO on 6 threads, alone on the machine): 71,188 nodes, 406,080 tetrahedra, certified (P4, P6); 17 of 17 points converged, one bias sub-step; build 5 s, sweep 892 s (FC10: 288,729 nodes, build 454 s, sweep 3,099 s). I_D 1.559e-10 A at 0 V to 1.575e-5 A at 0.8 V (FC10: 1.94e-10 to 1.37e-5 A, so the swept mesh is 15% higher in on-current and 20% lower in off-current); V_T (max g_m) 0.227 V; swing 64.3 mV/dec (FC10 about 62 mV/dec from 0 to 0.15 V; FC10 recorded no V_T). No reference exists for this shape, so neither run is known to be the closer; on the reduced fin (UC8) the swept mesh is within about 1% of the reference in on-current, the generated one 13.9% low.
- **UC11 met.** The gated fin (two straddled stations, a partial gate patch, P4 and P6), Release, idle machine,
  twice: 13,086 nodes (1,454 per section × 9 lines) in 0.90 s; **116,508 nodes** (5,548 × 21) in **8.0 s** (UD13's
  budget 10 s per 1e5).
- **A 2D generator defect found (not fixed; for review).** At (h, s) = (0.06, 0.04), (0.035, 0.02) and (0.03, 0.02)
  on the gated fin's cross-section (no layers), `generate_mesh(Geometry2D)` refuses its own result: "boundary facet
  is not Gabriel (P2)". The other 17 of 20 (h, s) pairs tried mesh. The same inputs fail identically on `fea05e5`
  (before Units 53 and 54), so the defect predates this unit; it refuses rather than returns a bad mesh. UC11 avoids
  those inputs (noted in the test). The cause is not diagnosed.
- **UC12, hand mutants (one per design row; Release, the mesh, deck and io tests):** UD1 a prism's third
  tetrahedron on the wrong diagonal (killed: UC3, UC4, the certificate); UD2 the lateral patches dropped from the 2D
  region signature (killed: the gated fin's patches); UD3 the region-change straddle not required (killed: the
  refusal); UD4 the void rule off (killed: the refusal); UD5 a parallel interface's upper straddle line at 0.6 s_i
  (killed: six layer cases); UD6 the mitre's kappa with 1 − cos (killed: the tapered fin; right angles unchanged);
  UD7 refinement allowed into a lattice (killed: seven cases); UD8 P6 skipped (killed: two); UD9 `max_nodes` against
  the section's nodes (**survived**: the test's cap of 50 was below either count; closed: the cap at the mesh's own
  count passes and one fewer is refused); UD10 the deck's layer `straddle=` ignored (killed: the deck tests); UD11
  the digest without a layer's straddle (killed: the digest test); the mesh file's corner blocks written as the
  layer count (**survived**: the round-tripped mesh has neither; closed: distinct counts saved and read back in
  place); UD12 U14's layers given the signed layers' exemption rule (**survived**: no test had a triangle of layer
  nodes outside every lattice; closed by a new UC1 pin, two U14 layers near a corner, 215 nodes, digest
  `1c2949ce70fe860a` measured identical on `fea05e5`; the mutant gives 210 nodes); UD13 the cost row has no output
  to mutate (argued equivalent; measured by UC11).
- Suites: Release 128/129 (618 s): `solve_mosfet` crashed (SEGFAULT, after 577 s, no assertion output); re-run alone it passes (12 passed, 1 skipped, 561 assertions). This slice changes no solver or library code (tests, docs and one deck only), and no such crash is recorded before; not reproduced, cause unknown, reported. Debug 118/118 (1,141 s).
- **CI (owner, 2026-10-31: "raise the Debug job's timeout-minutes from 60 to 90").** PR #74's last Debug job (run 41176925433) reached the 60-minute limit with 114 of 116 tests passed and none failed, and was cancelled; the previous run took 50 minutes. `.github/workflows/ci.yml` now gives the Debug preset 90 minutes and Release 60. Debug growth this slice, re-timed serially: `mesh.generate_swept` 32.2 s (+0.8 s, the `max_nodes` check), `deck.build` 15.1 s, `mesh.generate2d` 2.6 s; Unit 54 at about 56 s of UC13's 60 s.
