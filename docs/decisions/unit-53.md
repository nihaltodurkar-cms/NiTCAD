<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 53 — FinFET geometry and meshing (proposed 2026-10-29; on `mesh/finfet-geometry`)

**Owner request:** "finfet-geometry-scope implement" (2026-10-29): ARCHITECTURE.md 11.3, "FinFET geometry and meshing".
The branch starts from `b42253e` (the scope text, `docs/finfet-geometry-scope`, not yet on `main`), itself on `main`
`cb3f9de`; Unit 52's DECISIONS section (PR #73) and this one insert at the same place, a text-only merge.
**Owner rulings (2026-10-29, before the design):** (1) the deck describes the shape as a **sectioned extrusion**; the
general complex is API-only; (2) device gate: a **reduced FinFET** against a 2× refined tensor reference, Unit 16's
k = 1.5 rule; (3) **one boundary patch may span several facets** (one wrapped gate contact); (4) **item 4 is in this
unit**: patch ends and patch edges on interface traces are straddled, not refused. So the unit covers 11.3 items 1–4;
items 5 (3D boundary layers), 6 (doping junctions) and 7 (external import) stay planned.
The standing conditions hold: no change to any existing mesh, device or result; criteria fixed before code; Release
and Debug suites at every commit, `/W4 /WX`; hand mutants; nothing pushed or merged without authorization.

**Sources.**
- Shewchuk 1998 (SoCG, pp. 86–95), already read for U6: its input is a piecewise linear complex (vertices, segments,
  planar facets, not necessarily convex) and its output a *conforming* Delaunay tetrahedralization with radius–edge
  ratio ≤ 2 when no input angle is below 90°. To verify for this unit: that the hypothesis is stated for non-convex
  complexes with reflex (270°) dihedral angles, which the paper's figures suggest; until checked, "cited, to verify".
- Gärtner–Si–Fuhrmann 2010 (M3 b, M6): boundary-conforming Delaunay meshes for polyhedral domains; whether the
  admissibility of the clipped Voronoi dual (every volume > 0, boundary areas ≥ 0) is stated for non-convex domains is
  to verify. It does not matter for correctness: P3 is checked on every mesh.
- Legacy: none. The legacy unstructured 3D path (`unstructured_dd3d.py`) imports meshes; it has no generator.

**Evidence labels** as in Unit 16: proven, cited (to verify), heuristic (capped and checked after the fact), limitation.

**Design:**

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| FD1 | **The complex** (`mesh::PolyhedralGeometry`, API): vertices; facets (planar, simple polygons of vertex indices, at least three, no repeated vertex); regions, each a closed polyhedron given as a list of facets. A facet lies on one region (an outer boundary facet) or two (an internal facet). Rules, exact (`invalid_input`, index): every facet planar (orient3d) and simple (exact 2D tests in a projection); every region's facets a closed, oriented 2-manifold of non-zero volume; facets meet only along shared edges or at shared vertices (pairwise exact intersection tests; facet counts are hundreds); the domain connected. An internal facet between two regions of the same index is not an interface (dropped); between different indices it is an interface. | Exact predicates (M4). | A real FinFET is a union of non-convex regions; the complex is the standard input (TetGen's PLC). | O(F²) facet tests; fine at hundreds of facets. |
| FD2 | **Sectioned extrusion** (ruling 1; `mesh::ExtrudedGeometry`, deck and API): an axis (x, y or z); strictly increasing stations t₀ < … < t_k; a **base tiling** of the cross-section plane into named cells (simple polygons, a conforming tiling as Geometry2D's polygons); for each cell and section (t_i, t_{i+1}) a region index or *void*. It builds the complex with no polygon Booleans: lateral facets are cell edges × section; station facets are cells, an interface where the two sections' regions differ, internal where equal, boundary where one side is void. Adjacent coplanar facets of the same kind (same region pair, or same region and patch on the boundary) merge into one facet, so their shared edges are not features (FD7). A tapered fin is a trapezoid cell; the wrapped stack, spacers and raised source/drain are cells assigned per section. | Every vertex is a cell vertex at a station: exact coordinates, nothing is computed by intersection. | Covers the 11.3 shape (tapered fin, stack wrapped on three sides, spacers, raised S/D, STI/BOX) with exact inputs. | Faceted or rounded spacer *profiles along the axis* are not expressible; they need the API complex. |
| FD3 | **Boundary patches** (ruling 3): every outer facet carries one patch name; a patch may name any number of facets, and a node on an edge shared by two facets of one patch gets the sum of its two areas. M9's tiling of one facet into patches stays for partial facets. Patch names unique. | Under P2 a node's boundary area on a facet is its 2D Voronoi cell among the facet's nodes (M9's argument), so a multi-facet patch's area is the exact sum. | One wrapped gate contact instead of three. | A device consumer that assumed one face per patch: none found; checked in slice 2. |
| FD4 | **Meshing** (`generate_mesh(PolyhedralGeometry, MeshOptions3D)`): Shewchuk's conforming Delaunay refinement of the complex, extending `generate3d`. Every point inserted in one Delaunay tetrahedralization of a bounding box; subsegments are the complex's boundary edges, subfacets the triangles of each boundary facet's 2D conforming Delaunay triangulation in facet coordinates (non-convex facets: their edges as constrained segments, split in-plane when encroached). As now, a subsegment or subfacet is encroached when it is missing or a vertex lies strictly in its diametral ball, and encroached ones are split first (concentric shells at vertices), so missing boundary is recovered by Steiner points, never by constrained faces. Tetrahedra: radius–edge and size rule as now; a circumcentre that encroaches is not inserted, one outside the domain that encroaches nothing is exempt and counted. **Extraction:** flood fill from the bounding box's tetrahedra across triangles that are not subfacets; what is reached is outside. Then the flat-tetrahedron drop, and the checks: every boundary triangle on one facet, the volume and every facet's area equal to the complex's within 1e-12, P0–P3, P5. | Conforming, not constrained: P2 (every boundary triangle Gabriel) is what the box method needs (M3 b). Flood fill is exact once every subfacet is present. | Reuses Unit 16's tested refinement; no constrained tetrahedralization. | More Steiner points near small angles than a constrained mesher. |
| FD5 | **Termination:** input angles (between adjacent facets and segments, reflex included) of at least 90°: cited (Shewchuk 1998, FD sources); below 90° (a tapered fin's 85° foot, its 95° top is fine) heuristic, capped by `max_nodes` and the refinement iteration cap (`resource_exhausted`), as M11. | Cited / heuristic. | The theorem's hypothesis, not a refusal rule (M6). | A heuristic input may stop at the cap, with the cap named. |
| FD6 | **Region of a node:** exact point-in-region by ray parity against the region's facet triangulations (orient3d), with the ray's direction perturbed symbolically so it never meets an edge or vertex; boundary nodes take their facet's region. A node on an interface facet is a generator bug (the band forbids it) and fails with `degenerate_mesh`. Seeds are not needed: regions are geometry. | Exact. | Replaces the slab keys of planar interfaces. | O(nodes × facets of the region), measured (FC12). |
| FD7 | **Interfaces, bends and junction lines** (item 3): the interface pieces are the maximal coplanar unions of interface facets of one region pair. Features: **lines** (a piece's edges: a bend between two pieces of one pair, a junction of three or more regions, a *trace* where a piece meets the outer boundary) and **points** (where two feature lines meet, or a line meets a domain edge). Nodes (M7, U15 extended): on each piece, inset by 1.5 s from its feature lines, a pair set P from the 2D generator at spacing min(h, 4 s), boundary nodes graded 1.5 s, 3 s toward the inset polygon's corners, and for each p the nodes p ± (s/2) n; traces as now (pairs at ± s/2 within the face, graded from each end); no corner or bisector nodes anywhere (U15). Band b = √(R² + s²/4) around every piece: no other node, an encroached subfacet within b of a trace gets a new trace pair, tetrahedra of pair nodes only exempt (6d). Rules (`invalid_input`, the facet index): every dihedral sector at a feature line at least 60° (U12 extended; the tapered fin's 85°/95° pass); each piece's feature lines at least 6 s long; unrelated features at least 4 s apart; a trace only where the boundary is planar across it (two coplanar boundary facets), so an interface ending on a boundary edge is refused, as the 2D generator refuses one at a domain corner. | Proven away from feature zones (M7's mirror argument on each piece); a 90° bend is exact by mirror symmetry in 2D (U15) and expected so in 3D, measured; junction lines and points heuristic, judged by P4b. | The 2D rules that passed (U12, U15), lifted to lines; nothing new is assumed. | Interface geometry within 2 s of a feature is accurate to s only (U4). |
| FD8 | **Patches at traces** (item 4, ruling 4): on a boundary facet an interface trace is not a mesh edge but a line straddled by its trace pairs. A patch edge may (a) **lie on a trace**: the edge is not a subsegment, each face node belongs to the patch by its side, and by mirror symmetry the 2D Voronoi boundary between the pairs is the trace; (b) **cross a trace** at a point: a feature point on the face (sector at least 60°, corner zone 2 s, patch-edge nodes graded toward it as patch ends in 2D); (c) **span regions**: a patch over facets of different regions is the exact union. A patch edge elsewhere stays a constrained subsegment, exact (M9). The cut-face refusal of planar interfaces ("a face that carries patches", 6d) is lifted the same way. | Proven outside zones (M7 in 2D on the face, U2's argument); heuristic in zones, judged by P6. | Contacts on the fin's end faces and on raised source/drain beside the spacers. | A patch area is exact only to s × (zone length) near trace ends. |
| FD9 | **P4 extended** (`check_interfaces`): a declared piece may be a planar, simple, non-convex polygon. For a face F of a cross-region edge, deviation(F) is the largest distance of V(F) from the piece's plane when the 2D hull of V(F)'s projection lies inside the polygon (exact 2D tests): then the nearest point of the polygon to every point of F is on the plane, and plane distance is linear, so V(F) bounds it. Otherwise the piece's convex decomposition is used as now (min over pieces of max over V(F)). Corner zones only around feature lines and points, not around the decomposition's internal edges. Acceptance P4a and P4b unchanged; the certificate's form unchanged, its interfaces digest covering the new pieces. | Proven (convexity of the plane distance; the hull test is exact). | Without it every internal edge of a non-convex piece (the U-shaped spacer face at a gate edge) would become a corner zone, losing P4a there. | Two code paths in the check, both mutated. |
| FD10 | **P6, patch boundaries** (new check, in the certificate): on each boundary facet, every 2D dual edge between face nodes of different patch membership (the boundary piece of the dual face of a boundary edge) lies on its declared patch boundary within P4a's tolerance outside the zones and within s inside them; every patch's area against its geometric area reported. `Device::create` requires it wherever FD8 straddles a patch edge (always a multi-region device, which already needs a certificate). | Check, not claim (M12's practice). | FD8's patch edges are not mesh edges, so their exactness must be verified like interfaces. | — |
| FD11 | **Options:** `MeshOptions3D` unchanged (h, per-region h, size boxes, s, radius–edge, `max_nodes`). | — | — | No anisotropic layers (item 5): FC9's risk (below). |
| FD12 | **Existing paths:** `generate_mesh(Geometry3D, …)` (convex, planar interfaces) stays and may share code with FD4, but every mesh it generates is bit-identical to `main` (FC1). The 2D generator and the tensor path are untouched. | — | No existing deck, mesh file or result changes. | Shared internals must keep the convex path's insertion order. |
| FD13 | **Deck** (ruling 1): `[mesh] type = extrude3d`; `axis = x|y|z`; `stations = <lengths>`; `cell <name> = <polygon in the cross-section's two coordinates>`; meshing keys as `generate3d` (`max_size`, `interface_resolution`, `box`, `radius_edge`, `max_nodes`). In a `[region]`: `extrude = <cell names> <t_lo>:<t_hi>` (t at stations), several lines allowed; a (cell, section) given twice is an error; unassigned ones are void. In a `[contact]`: `patch = <name>` as now; patches are named in `[mesh]` by `patch <name> = box <lo> <hi>`, the whole boundary facets inside the closed box (a box selecting none is an error; facets are never clipped). In an extrusion every trace and every patch edge is a facet edge, so FD8's cases (a)–(c) all arise from whole-facet patches. | Exact selection: nothing is computed from the box except containment. | Users write a FinFET as cells and stations, as its cross-section drawings are drawn. | A patch through a facet's interior needs a cell edge or station there. |
| FD14 | **io, run record, app:** generated-mesh files as M14 (simplices, regions, patches, certificate with its P6 report, a generation digest over the complex and options); the run record names the mesh type and its digest; `nitcad check` reports node, tetrahedron, sliver and exempt counts; the application renders it as any generated 3D mesh. | — | Reloaded meshes stay renderable and certified. | — |
| FD15 | **Determinism and cost:** BRIO with the fixed seed as M13; the same complex and options give a bit-identical mesh in one build. Budget, measured not promised: 1e5 nodes ≤ 30 s Release (twice M13's convex budget, for boundary recovery and node location). ctest meshes ≤ 2e3 nodes as M13, larger `[.benchmark]`. | — | — | — |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| FC1 | Nothing changes: every existing test and pin passes; every mesh `generate_mesh(Geometry3D)`, the 2D generator and the tensor path produce in the suites has the same graph digest, region digest and certificate as on `main`; every existing deck's run identity and device file byte-identical. |
| FC2 | Refusals: every FD1, FD2, FD7, FD8 and FD13 rule refused with its code, its index or deck line, and its message (one test each). |
| FC3 | Non-convex single region: an L-prism, a U-channel, a fin on a slab (T), a slab with a notch (reflex edges), a tapered fin (85°), a facet whose boundary is non-convex: P0–P3 and P5 hold; every boundary triangle Gabriel and on one facet; the volume and every facet's area equal the complex's within 1e-12; the radius–edge ratio ≤ the bound outside exempt tetrahedra; two runs bit-identical. |
| FC4 | Multi-facet patches: a patch over the three facets of a wrapped surface has, at every node, the sum of the per-facet 2D duals' areas within 1e-12, and its total equals the three facets' geometric area within 1e-12. |
| FC5 | Linear field on non-convex domains: for u = a·x_k (k = each axis on the L, U and T shapes; the fin's axis on the tapered fin), Dirichlet on every facet normal to x_k and Neumann elsewhere, conjugate gradients on the generated dual give u at every node within 1e-12 of its range. |
| FC6 | Interfaces (P4 certified, P4b ≤ s): a fin with a liner wrapped on three sides (90° bend lines), the same tapered (85°/95°), a liner foot on the slab (three-region junction line), a junction point where those lines meet a gate-edge plane, a trace on a flat boundary, a 60.5° sector. Every flat deviation within P4a; every corner deviation measured, recorded and pinned with a regression bound taken from the measurement. Refused: a 59.5° sector, an interface ending on a boundary edge, features within 4 s. |
| FC7 | Patches at traces (FD8, P6): a patch ending on a trace, a patch edge crossing a trace, a patch spanning two regions: certified (P6), the patch boundary's flat deviation within P4a, its corner deviation ≤ s and pinned, each patch's area error reported and ≤ s × the zone length. |
| FC8 | Cross-section convergence (the gate for FD4 + FD7 together): a gated fin uniform along its axis (tapered fin, liner and high-k wrapped, the wrapped gate as one patch, ohmic contacts on the fin's end faces, zero drain bias, gate at 0.5 V), 3D generated meshes at h, h/2, h/4 against a 2D `generate2d` equilibrium solve of the cross-section 4× finer than the finest: discrete H1 potential order ≥ 0.9, L2 ≥ 1.5, the gate charge per unit length order ≥ 0.9 (Unit 16's thresholds); every mesh P4 certified. |
| FC9 | Reduced FinFET (ruling 2): the tensor deck's shape (rectangular fin, box domain, spacer fill, one gate patch over three faces) with a coarser tensor mesh whose 2× refined tensor reference solves within about an hour (PARDISO, Release, `[.benchmark]`). The generated mesh has at most the reduced tensor mesh's nodes. For I_D at every V_G point above the current resolution, V_T, max g_m and the swing: \|Q_gen − Q_ref\| ≤ 1.5 \|Q_tensor − Q_ref\| + f, f the reference's Richardson estimate. Every point converges (step cutting recorded if used). **Risk stated before the run:** the 2D MOSFET met this rule only with a structured layer (U14), which is item 5, outside this unit; if FC9 fails it stays in the suite as `[!shouldfail]` with its measurement, for an owner ruling (pull item 5 in, or change the gate), never a loosened threshold. |
| FC10 | The shaped FinFET deck (`decks/finfet_3d_shaped.deck`: tapered fin, wrapped stack, spacers, raised source and drain, STI/BOX, one gate contact, source and drain contacts ending on traces): meshes, is certified (P4, P6), and its Id–Vg sweep converges at every point; nodes, slivers, mesh time and solve time reported; no accuracy gate (no reference). |
| FC11 | Deck, io, run record, app: the extrusion vocabulary checked by `nitcad check` with line numbers; a generated mesh round-trips with its patches, regions, certificate and P6 report; the run identity covers the complex; `nitcad plot` and the application show the mesh and fields. |
| FC12 | Cost reported: mesh time and node location time at about 1e4 and 1e5 nodes, against FD15's budget. |
| FC13 | Hand mutants: at least one per design row (FD1–FD15), each survivor closed by a new test or argued equivalent. |
| FC14 | Release and Debug suites pass at every commit, `/W4 /WX`; the Debug suite grows by ≤ 120 s or cases move to `[.benchmark]` (M13). |

**Slices** (each a local commit with both suites green):
1. The complex and the extrusion builder, validation and node location (FD1, FD2, FD6; FC2 part).
2. Non-convex single-region meshing and multi-facet patches (FD3–FD5, FD12; FC1, FC3, FC4, FC5).
3. Interfaces with bends, junction lines and points; P4 for non-convex pieces (FD7, FD9; FC6).
4. Patches at traces and P6 (FD8, FD10; FC7).
5. Deck, io, run record, application, the two FinFET decks (FD13, FD14; FC11).
6. Gates and cost: FC8, FC9, FC10, FC12, mutants FC13, ARCHITECTURE 6.11 and 11.3.

**Owner questions (FQ):**
- FQ1: approve FD1–FD15 and FC1–FC14 as written?
- FQ2: FC9's pre-stated failure path (`[!shouldfail]` and a ruling, no threshold change): agreed?
- FQ3: FC8's bias point (0.5 V gate, inversion) only, or a C–V sweep (−0.5 to 1 V) with the order gates at each point?

**Owner answers (2026-10-29), design and criteria approved:** FQ1 approved as written; FQ2 a failed FC9 stays
`[!shouldfail]` with its measurement for a ruling, the threshold never loosened; FQ3 FC8 at the 0.5 V gate only.

**Slice 1 (2026-10-29): the complex, the extrusion and node location** (`include/NiTCAD/mesh/complex.hpp`,
`src/mesh/complex_index.hpp`, `src/mesh/complex.cpp`; `tests/mesh/complex_test.cpp`, 5 cases):
- `check_complex` and `detail::ComplexIndex` implement FD1 exactly (orient2d/orient3d only). Facets are triangulated by
  ear clipping in a one-to-one coordinate projection, a strictly convex ear holding no other vertex, so collinear facet
  vertices never give a flat triangle. The pairwise facet test is the triangle-pair rule "the intersection lies in the
  convex hull of the shared vertices", checked edge against triangle both ways; the argument (the intersection of two
  triangles is convex and its boundary lies on their edges) is in `complex.cpp`.
- **Added to FD1, for review:** volumes with disjoint surfaces must not nest (a vertex of one, not shared, located
  against the other). FD1's facet test alone does not see a volume inside another.
- FD6: `locate` gives on a facet (closed), inside a volume, or outside. A segment that meets a triangle's edge or vertex
  is degenerate and is replaced by the next of a fixed splitmix64 sequence of far points, so the result is exact and
  deterministic. Checked on a quarter-step lattice through two cubes, which hits every tie with their planes and edges.
- FD2: `extrude` checks the tiling itself (whole shared edges in opposite directions, no vertex inside another cell or
  inside another cell's edge, no crossing, no collinear overlap), so its messages name the cell. It then builds one
  volume per non-void (cell, section) and passes the result through `check_complex`. A patch box takes the boundary
  facets inside it that no earlier box took. A box that takes none is refused.
- Not testable by construction, kept as defensive rules (mutated in slice 6): "not orientable" (an embedded closed
  edge-manifold surface is orientable) and "zero volume" (a flat closed surface fails the facet test first).
  "Pinched at a vertex" is reachable only through a connected surface that touches itself at a vertex; no fixture
  builds one yet.
- Suites: Debug 109/109; Release 119/120 in the full `-j 6` run. `app.live_check` failed once; its output was cut off
  by my log filter, so the failing case is unknown. It then passed alone, and 15 times in a row with the app tests
  under `-j 6`. It does not link the new code (mesh layer only). Unit 50 fixed one race in this test (XG1), so this
  may be another; it is recorded open, not retried away.

**Slice 2 (2026-10-30): non-convex single-region meshing and multi-facet patches** (`generate_mesh(PolyhedralGeometry,
MeshOptions3D)` in `src/mesh/generate_complex.cpp`; `tests/mesh/generate_complex_test.cpp`, 5 cases):
- FD12 as built: a separate generator file. `generate3d.cpp` is not touched, so FC1 holds for the convex path by
  construction. The new file repeats the convex generator's refinement loop, about 400 lines. Merging the two is
  deferred until the interface code (slice 3) shows what they share.
- Faces are the maximal connected sets of coplanar boundary facets: non-convex and with holes (a fin's footprint on a
  slab top). Domain edges are chains of collinear facet edges where two faces meet. Patch edges are facet edges inside
  a face between different patches. A facet edge between facets of one patch in one face is not a mesh constraint.
- **Added to FD4, for review:** a domain-edge subsegment is also encroached when it is missing from the 2D
  triangulation of either face it bounds (a reflex edge of a non-convex face can be absent there while present in
  3D). Only tetrahedra whose centroid is in the domain are refined (exact `locate`). Both are needed only for
  non-convex input.
- **Added rule:** the boundary must be a manifold surface (an edge in four boundary facets is refused, `invalid_input`).
- Extraction: a flood fill from the bounding tetrahedra; it stops at triangles inside one face and skips tetrahedra flat
  on a face. Then: every boundary triangle inside one face; the volume (1e-12); every face's area (1e-12, new); then
  `make_dual_mesh`.
- Measured (Release, `max_size` 0.45–0.8 on O(1) shapes): L-prism 346 nodes, U-channel 725, pocket 886. FC5 linear-field
  deviation at most 1.6e-15 relative on every axis of the L, U, T and pocket shapes, and on the tapered fin's axis.
  FC4: the three-facet gate patch's node areas match the per-plane 2D duals within 1e-12 (8.0 total). FC3: P0–P3 and P5
  (`make_dual_mesh`), plus a Gabriel recheck in the test, on all six shapes, including the 86° tapered fin (heuristic)
  and a non-convex L facet given directly.
- Suites: Release 121/121. Debug 109/110: `app_window.self_test` failed on the WG1 gate (the message loop's longest
  gap 110.7 ms, limit 100 ms, Microsoft Basic Render Driver). Run alone it failed once more, then passed twice in a
  row. This is the open WG1 spike (ARCHITECTURE 11.3, cause unknown); the application does not call the new
  generator.

**Slice 3 (2026-10-30): interfaces with bends, junction lines and points; P4 for non-convex pieces**
(`src/mesh/complex_interfaces.{hpp,cpp}` (FD7 layout and rules), `src/mesh/polygon2d.{hpp,cpp}` (the 2D polygon
helpers, moved from `complex.cpp`), `generate_complex.cpp`, `interfaces.cpp` (FD9), new exact predicates
`dot_sign_3d`, `dihedral_dot_sign`, `angle_below_60_3d`, `dihedral_below_60`;
`tests/mesh/generate_complex_interfaces_test.cpp`, 5 cases, and one FD9 case in `interfaces_test.cpp`):
- FD7 as built. A piece is a union-find group of coplanar interface facets of one region pair; its loop is the edges
  used once. Feature lines are classified as trace or interior. Every dihedral sector about a line is sorted exactly
  (orient3d and the sign of the dihedral dot product) and checked against 60° exactly; the exterior sector of a trace
  is skipped. Feature points: crossings of domain edges (a single piece, both arms on the edge's two faces, 60° to the
  edge, no domain-edge vertex within 4 s) and trace junctions inside a face (60° in the face).
- Placement as built:
  - each piece's pair set is the 2D generator's nodes in its 1.5 s inset at min(h of both regions, 4 s);
  - each trace gets pairs at ± s/2 in its face;
  - each crossing gets a pair at ± s/2 along the domain edge, and the domain-edge segment across the crossing is
    locked (never split, never encroached).
  An encroached subfacet whose circumcentre is in the band gets a new trace pair at the projection on the nearest
  trace of its face. A tetrahedron of interface nodes only, or with its circumcentre in the band, is exempt. Node
  regions come from exact `locate` (FD6); a node on an interface facet is `degenerate_mesh`. Sizes per region: a
  point on an interface takes the smaller h.
- **Changed from FD7, for review:** FD7's "graded 1.5 s, 3 s toward the inset polygon's corners" is replaced by a
  **line lattice**. Rows along every feature line sit at multiples of 1.5 s in a sign-canonical coordinate along the
  line, and the trace pairs and the pair sets' inset-edge nodes use that lattice. So the pieces meeting at a line
  (a bend, the three of a junction) have aligned rows, and near the line the nodes are an extrusion of the 2D
  cross-section. The graded rule left the rows misaligned at the gated fin's T-junction (P4b failed); densifying the
  inset boundary made it worse. Crossings get no corner nodes (U15).
- **Added rules, for review** (cautious defaults; the owner was asked and has not yet ruled):
  - the 4 s separation also applies between pieces that do not touch (a layer thinner than 4 s; parallel pieces) and
    between a piece and a boundary face it does not meet;
  - refused: holes in a piece, a pinched piece, a piece too narrow for its 1.5 s inset, a trace at a domain corner, a
    feature line on a domain edge other than a single crossing;
  - refused for now: a patch edge meeting an interface (slice 4 lifts it, FD8);
  - an oblique non-convex piece always takes FD9's decomposition (its orthogonal projection is not exact).
- **FD9 as built, with one addition for review.** A polygon the convex rule accepts stays on the original path, byte
  for byte. Otherwise the polygon must be planar and simple in the projection dropping the normal's largest component,
  and it is ear-clipped. For an axis-aligned piece, the exact hull test (monotone chain, every hull vertex in the
  closed polygon, no polygon edge meeting the open hull) gives the plane distance. **Added:** a second hull test
  against the polygon grown outward by 1e-6 s (mitred) gives the plane distance plus the grown reach. Without it,
  the gated fin failed P4b: face vertices rounded 1e-16 outside a feature line fell to the decomposition (0.11). Such
  faces are in corner zones and judged against s, so the reach (≈ 1e-6 s) is immaterial. Otherwise the
  decomposition's min over parts. Direct test (`interfaces_test.cpp`): an L-shaped piece and its notch on a lattice
  whose Voronoi lines fall on the L's edges is accepted; without the notch piece it is refused; a bowtie is refused.
  Two hand mutants checked now (FC13 proper is slice 6):
  - hull path disabled: the L case fails, and so does one `[complex_interfaces]` assertion;
  - hull test always true: the notch refusal fails.
- **FC1 by construction, not by a digest dump.** `generate2d.cpp`, `generate3d.cpp` and the tensor path are not
  touched; `predicates.cpp` only gains functions. In `interfaces.cpp`, every polygon `main` accepts runs the same
  operations in the same order (`parts` empty), so every existing report and certificate is unchanged. Only the
  refusal text changed, "convex" → "simple", and no existing test matched it.
- **Fixed in slice 2's code:** with interfaces every node is located, and a node on an oblique boundary face can
  round just outside the domain. It now takes the region beside the nearest boundary facet within 1e-9 of the
  domain's extent. Meshes without interfaces keep slice 2's single region and never call `locate` (the first full
  Release run failed two slice-2 cases on this, the tapered fin and FC5).
- FC6 measured (Release, s = 0.1, liner t = 0.7, h = 0.5). Pins are regression bounds from these measurements:

  | Mesh | Nodes | Flat deviation (P4a ratio) | Corner deviation | Pin | Slivers |
  |---|---|---|---|---|---|
  | Liner fin, 90° bends, junction lines at the feet | 1469 | 0 | 0.286 s | ≤ 0.32 s | 71 |
  | Tapered fin, 85°/95° | 1456 | 1.4e-3 | 0.141 s | ≤ 0.16 s | 115 |
  | Gated fin, inverted-U pieces, junction points | 3417 | 0 | 0.411 s | ≤ 0.45 s | 262 |
  | 60.5° sector | 725 | 2.2e-3 | 0.446 s | ≤ 0.49 s | 50 |

  For the liner fin, 1032 of its nodes are interface nodes. The worst corner faces: on the liner fin, at the end
  face where the liner's outer corner line meets it; on the gated fin, at the junction point of a liner corner line
  and the gate-edge plane. Both runs of the liner fin are bit-identical.
- The gated fin is 3417 nodes, over FD15's 2e3 ctest guideline. It stays in ctest because it is the only FC6 junction
  point fixture; its Debug time is below, under FC14's budget.
- FC2/FC6 refusals tested: no resolution, 59.5°, an interface ending on a domain edge, two interfaces ending along one
  line, a thin layer, near a parallel face, a short trace, a short feature line, a hole, unrelated features within
  4 s, a feature within 4 s of a domain edge, a domain-edge vertex within 4 s of a crossing, a patch edge on a trace.
  **Not yet reached by a fixture, kept as defensive rules (mutated in slice 6):**
  - a pinched piece;
  - a trace at a domain corner;
  - a feature line on a domain edge (more than a single crossing);
  - traces meeting in a face under 60°;
  - a trace meeting a domain edge under 60°;
  - a piece too narrow for its inset.
  An extrusion's cell grid reaches these awkwardly, and some may need the API complex.
- Suites:
  - Release 122/122. The first full run failed the two slice-2 cases above; after the fix, all pass.
  - Debug 111/111, total 912 s (4 jobs).
  - FC14: the new `mesh.generate_complex_interfaces` takes 85.0 s in Debug, and `mesh.interfaces` 17.8 s including
    the FD9 case (0.56 s in Release). That is within the 120 s growth budget but uses most of it, mainly the gated
    fin. Slice 4's fixtures must move a case to `[.benchmark]` if they would exceed it.

**Owner ruling on slice 3 (2026-10-30), "choose recommended":** every slice-3 item marked "for review" is accepted as
built:
- the line lattice replaces FD7's graded rows;
- the added rules stand: the 4 s separation for thin layers, parallel pieces and nearby boundary faces; refusals of
  holes, pinched and too-narrow pieces, traces at domain corners, and feature lines on domain edges other than a single
  crossing;
- oblique non-convex pieces always take FD9's decomposition;
- FD9's grown-polygon path stays.

**Slice 4 plan, owner ruling (2026-10-30), "choose recommended" (SQ1–SQ6):**
- SQ1: `generate3d.cpp` keeps its refusal "an interface cuts a face that carries patches"; FD8's sentence lifting it
  is amended to the complex generator only (FD12, FC1: the convex generator is not touched). The same geometry is
  meshed through the complex API.
- SQ2: `SimplexMesh` gains an optional per-corner patch (`facet_corner_patch`, D entries per listed facet; empty means
  `facet_patch` for every corner). `make_dual_mesh` gives each corner's boundary piece to its corner's patch. Under
  FD8 (a) a boundary triangle straddles the trace, and the boundary between a mirror pair's pieces is their
  bisector, the trace, so the split is exact; P6 judges the rest. Empty, the importer runs as before (FC1).
- SQ3: FD10's requirement is enforced when the certificate is issued: `certify_interfaces` refuses simplices with
  per-corner patches unless a patch-boundary set is given, so every certificate of such a mesh carries P6.
  `graph_digest` covers the boundary patches, so `Device::create`'s digest check keeps the certificate tied to them.
- SQ4: until slice 5 (FD14), writing a generated mesh with per-corner patches or a P6 report is refused, never
  written without them.
- SQ5: FD8 (b) covers every point where a constrained patch edge meets a trace: a crossing, and a patch edge ending
  on a trace (a T, the straddled part continuing). Rules as a domain-edge crossing: sectors at least 60°, 4 s
  separation, not at a domain corner.
- SQ6: FC7's "zone length" is the length of declared patch boundary inside the 2 s zones.

**Slice 4 (2026-10-30): patches at traces and P6** (FD8, FD10; FC7, FC2 for the FD8 rules).
- **Baseline first:** the main merge `8f30ac0`, never built before its push, passes Release 124/124 (779 s) and Debug
  113/113 (1181 s, 4 jobs).
- **FD8 (a), as built.** In `generate_complex.cpp` a facet edge between two patches of one face that is also an edge
  of an interface facet is recorded as a straddled patch edge and left out of the chains, so it is never a
  subsegment. The trace pairs and the line lattice are unchanged. On a face with a straddled edge each boundary
  triangle's corner takes the patch of the face triangles (closed, in face coordinates) holding its node; when they
  disagree (a node on a constrained patch edge), the triangle's own patch (the centroid's, as before). The corners go
  into `SimplexMesh::facet_corner_patch` only when some corner differs from its triangle, so every mesh without a
  straddled edge is built and imported exactly as before (FC1 by construction).
- **FD8 (b), as built (SQ5).** In `read_interfaces` a patch chain through an interface loop vertex is a crossing
  (`InterfaceCrossing`, its position 0 or last for a patch edge ending on the trace). Refused: a domain chain or a
  second patch chain through the point; anything other than one piece whose two trace arms lie in the patch edge's
  face; a sector under 60 degrees between a trace arm and the patch edge; a patch-edge vertex within 4 s of the point.
  The generator puts nodes at +-s/2 along the patch edge (one at an end, inserted in the patch edge's one face), the
  segment across never split, as at a domain-edge crossing; on the trace through the point no lattice pair within
  0.75 s, and the point's position is kept among the trace's positions so refinement adds none within s / 2.
- **Added rule, for review:** a patch edge must be at least 4 s from every feature line it does not meet (the domain
  edges' rule extended to patch edges; their subsegments would otherwise be encroached by the pairs).
- **FD10, P6, as built** (`check_patch_boundaries`, `interfaces.hpp`): per listed boundary triangle, each edge whose
  corners carry different patches; the piece from the edge's midpoint to the triangle's circumcentre; deviation min
  over the declared boundaries of those two patches of max over the piece's two ends; zones of 2 s about declared
  centres; P4a's tolerance outside, s inside; patches meeting with no declared boundary refused. Reported: each
  patch's area error and its zone length (sampled at 64 points per zone radius; a report, not judged). The generator
  declares every patch edge (straddled and constrained) and puts zone centres at the ends of straddled runs, where
  they meet a chain, and where patch edges meet traces. It runs P6 only where a patch edge meets an interface
  (straddled or crossing), so every slice-1 to slice-3 certificate is unchanged.
- **SQ3 as built:** `certify_interfaces(…, const PatchBoundarySet* = nullptr)` refuses simplices with per-corner
  patches without a set; the set enters `interfaces_digest` only when given (existing digests unchanged).
- **SQ4 as built:** `save_generated_mesh` refuses per-corner patches or a P6 report; `save_device` refuses a
  certificate with a P6 report; both before anything is written (tested).
- **SQ1:** `generate3d.cpp` not touched.
- **FC7 measured** (Release, s = 0.1, h = 0.5; slabs 2 × 2 split at z = 1). Pins are regression bounds:

  | Mesh | Nodes | P6 pieces (in zones) | P6 flat / corner | P4 corner | Area error | s × zone length |
  |---|---|---|---|---|---|---|
  | Patch edge on a trace (a) | 350 | 28 (13) | 0 / 0 | 0 | 0 | 0.04 |
  | Patch edge crossing a trace, patch over two regions (b, c) | 536 | 0 | 0 / 0 | 0 | 0 | 0.04 |
  | Patch edge ending on a trace (a, b) | 539 | 26 (11) | 0 / 0 | 0.25 s | −0.00125 | 0.06 |

  The straddled boundaries are exact even in their zones: the trace pairs are mirrored to the trace's ends. In the
  T case the P4 corner deviation and the area error are likely the patch edge's end node, s / 2 from the trace with no
  mirror (not traced). In case (b) the patch edge stays constrained, so no triangle straddles a boundary and P6 has
  no pieces; its report (areas, zone length) is still issued.
- **Refusals tested (FC2):** under 60 degrees (a cell edge at 53 degrees on the end face), a patch-edge vertex 0.2
  from the meeting point, a patch edge 0.2 from a parallel trace. The first fixture used 0.3 and 1.7, whose chain
  was not exactly collinear and split at the crossing (another refusal fired); exact coordinates fixed it.
  Not reached by a fixture, kept as defensive rules (mutated in slice 6): a domain edge or second patch edge through
  the point; trace arms in another face.
- **Direct tests:** `simplex_mesh_test.cpp` (corner patches equal to the facet's give the same graph digest; a corner
  patch moves exactly its corners' pieces; size and index refusals); `interfaces_test.cpp` (P6 on a mirrored
  straddled trace: exact, areas exact, zone length 4 s; a face node moved 0.002 refused outside a zone and accepted in
  one; no declared boundary refused; certification without a set refused; the set changes the digest).
- **Hand mutants now** (FC13 proper is slice 6), each killed:
  - corner patch always the triangle's: 2 `[complex_patches]` cases fail;
  - P6 zones never matched: 2 test cases fail (of the `[complex_patches]` and direct P6 cases run);
  - certification without a set allowed: 1 test case fails.
  A first form of the second (zone radius 0) hung the zone-length sampler (an infinite step count), not a valid
  mutant; replaced.
- Suites (2026-10-30, after the mutants were restored):
  - Release 125/125, 723 s (4 jobs); `mesh.generate_complex_patches` 0.50 s.
  - Debug 114/114, 995 s (4 jobs; the baseline took 1181 s, so totals vary run to run).
  - FC14: the new `mesh.generate_complex_patches` takes 9.0 s in Debug, and `mesh.interfaces` grew from 17.8 s to
    26.1 s (the direct P6 case, brute-force Delaunay). Slice 4 adds about 17 s, so Unit 53's growth is about 102 s of
    the 120 s budget (slice 3's 85 s plus this). Slices 5 and 6 have about 18 s left: their larger fixtures go to
    `[.benchmark]`.

**Slice 5 plan, owner ruling (2026-10-30), "choose recommended" (SQa–SQd):**
- SQa: one mesh and device file format; the corner-patch array and the P6 report are new optional records, so files
  without them read and write byte-identically (FC1). The SQ4 refusals are lifted.
- SQb: `patch <name> = box <lo> <hi>` selects the whole boundary facets inside the closed box (FD13); a box that
  selects none is an error with its deck line.
- SQc: the shaped deck keeps `finfet_3d.deck`'s dimensions and materials, the fin tapered about 5°, and its 50 mV
  Id–Vg sweep (FC10 has no accuracy gate).
- SQd: slice 5 builds and certifies both decks (the shaped one and FC9's reduced one) without solving; the solves are
  slice 6's FC9 and FC10.

**Slice 5 (2026-10-30): deck, io, run record, application, the two FinFET decks** (FD13, FD14; FC11).
- **Deck** (`src/deck/build.cpp`, `include/NiTCAD/deck/build.hpp`): `type = extrude3d` with `axis`, `stations`,
  `cell <name> = <u v>, ...`, `patch <name> = <x0 y0 z0 x1 y1 z1>`, the `generate3d` meshing keys (now one shared
  function), and `extrude = <cells> <t_lo>:<t_hi>` in a region. Refused with the deck line: a bad axis, no axis,
  stations or cell, an unknown cell, a range whose ends are not stations or not increasing, a (cell, section)
  assigned twice, `box`/`polygon`/`seed` in an extrusion's region, `extrude` in any other mesh type. A patch box
  that selects nothing is `extrude`'s own error (SQb), located at the `[mesh]` section (the library gives the box's
  index, not its line); the deck-level check written first was never reached and was removed.
- **io** (`src/io/files.cpp`, `device_file.cpp`): optional records, so older files are unchanged (SQa). The mesh file
  gains `facet_corner_patch` (key 15); the certificate record gains the P6 report (keys 4–7: counts, deviations,
  per-patch area error and zone length). The SQ4 refusals are gone. `load_generated_mesh` rebuilds the dual from the
  stored corner patches before comparing digests. Tests: a complex-generated mesh with corner patches round-trips
  (`generated_test`), and the device file round-trips a certificate with a P6 report byte for byte
  (`device_file_test`).
- **Run record:** not changed. The run identity already hashes the whole built mesh (points, edges, boundary patches
  and areas), so any change to the complex or its options that changes the mesh changes the identity (tested: a
  moved station, a renamed patch). FD14's separate "generation digest over the complex and options" is not built:
  the same mesh gives the same solve, and the digest would only add the options that do not change the mesh.
  **For review.** `nitcad check` now also prints the generated mesh's counts and the certificate (tetrahedra,
  slivers, exempt, smallest dihedral, largest radius–edge; P4 faces and P6 pieces), and `fail` prints the error's
  index and value.
- **Application:** the vocabulary (`deck_vocabulary.cpp`) and the section menu (`deck_sections.cpp`) know the new
  type and keys (the vocabulary test builds every key against an extrude3d fixture and checks the reference header
  in both directions); the workspace renders the mesh as any generated 3D mesh, unchanged.
- **A slice 4 defect found by the first real geometry, fixed:** a corner of a straddled patch boundary (the patch
  turns a corner where two interface traces cross on a face) is a P6 zone; it was not, because only run ends and
  chain meetings were. A single gate patch on a side face failed P6 by s / 2 at the corners. Fixed in
  `generate_complex.cpp` (a degree-2 vertex whose two straddled edges are not collinear is a zone centre); new test
  `complex patches: a patch corner where two traces cross`. P4b and P6 messages now say where (the edge's midpoint
  and the two patches).
- **Limitation found, not fixed:** two cells of one region that meet on a boundary face leave a complex vertex on an
  interface trace (the seam's end). The straddled construction assumes no node on a trace, and P6 then fails by
  s / 4 beside the vertex (reproduced with a slab split into two cells at a trace; one cell passes). The decks use
  one U-shaped cell per layer for this reason. A mesher that removed such vertices, or an `extrude` that merged
  same-region neighbours, would lift it; neither is built.
- **The decks, and a sizing finding.** With the reference deck's layers (0.5 nm SiO2, 1.5 nm HfO2) the generator
  needs s ≤ 0.083 nm: the interfacial layer's end ring at each station has 0.5 nm edges, and FD7 wants every feature
  line at least 6 s. That did not mesh within 200000 nodes (50 s to reach the cap). Both decks therefore use 1 nm
  SiO2 and 2 nm HfO2 (EOT about 1.2 nm) at s = 0.16 nm, which is **a change from SQc** (the dimensions of
  `finfet_3d.deck`), flagged for review.
  - `decks/finfet_3d_reduced.deck` (rectangular fin, box domain, one gate patch over three faces): meshes and is
    certified (P4 and P6) in 4 min 8 s Release: 129187 nodes, 753953 tetrahedra, 33707 slivers (dihedral below 5°),
    smallest dihedral 0.00°, largest radius–edge 41.0, P4 55768 faces (12611 in zones), P6 6244 pieces (150 in
    zones). Far over FD15's 30 s at 1e5 nodes and with poor tetrahedra: reported, not hidden.
  - `decks/finfet_3d_shaped.deck` (fin tapered 3.8°, the interfacial layer following it): **does not certify.** P4b
    fails at (30, 50.76, 11.775) nm, the top corner of the interfacial ring on the x = 30 plane: deviation 0.241 nm
    = 1.5 s against s = 0.16 nm. At s = 0.1 nm (21 min, the ring now 10 s thick) it fails at the same corner,
    (30, 50.85, 11.86) nm, by 0.151 nm = 1.5 s again, so the error scales with s and is not an overlap of zones:
    a corner defect of the tapered ring (right-angle rings, as in the reduced deck, pass). Not diagnosed further;
    the deck stays as the measurement, in a hidden `[!shouldfail]` test, for an owner ruling (FC10 is slice 6).
    Two layouts were refused first for FD7 reasons: an HfO2 outer face that bends where the buried oxide's
    interface ends (an interface may not end on a domain edge), and a void fill beside the gate (the interface
    plane then abuts a coplanar boundary face).
- **Suites** (2026-10-30, `/W4 /WX`, 4 jobs): Release 125/125 (739 s), Debug 114/114 (1119 s; slice 4's was 995 s, and
  totals vary run to run, so FC14's 120 s is not re-measured per test). Slice 5's new cases are small (the fixtures
  mesh in under a second); the two FinFET decks are hidden `[.benchmark]` cases.

**Owner ruling on slice 5's open items (2026-10-30), "fix them":** the shaped deck's P4b failure, the thicker
layers and FD14's generation digest are each fixed in code, not accepted as limitations.
- **FD14, the generation digest** (`6ec2286`): `mesh::generation_digest` (FNV-1a over the complex's vertices,
  facets with patch names, volumes with regions, and every `MeshOptions3D` field), stored in `GeneratedMesh` by
  the complex generator only (0 for the others), written as optional mesh-file key 16 (other generated files
  unchanged). The run record gains the mesh type and that digest as an optional record (key 5), set by the deck
  runner for every analysis on a complex-generated mesh; `nitcad check` prints it, the application's run panel
  shows it. It is not part of the input identity: `max_nodes`, for one, changes the digest and not the mesh
  (tested).
- **P4b at the tapered ring: a check bound, not a mesh defect.** Reproduced in 2.6 s by a 3-section miniature
  deck; the regression test `complex interfaces: a tapered liner ring at a gate edge` is a two-section reduction
  of it (0.7 s Release, 13.7 s Debug: FC14's Debug growth, about 102 s after slice 4, is to be re-measured in
  slice 6). At the station plane the interfacial ring (a
  non-convex piece) meets the coplanar HfO2 piece at oblique corners (about 94 and 266 degrees); a dual face
  there has vertices up to 0.04 s outside the ring's outer edge, so the exact hull test and the 1e-6 s grown
  polygon both fail, and the decomposition bound (min over ear triangles of max over V) is 1.1 s because the
  face straddles a diagonal of the decomposition, while the outermost vertex is 0.006 nm (0.04 s) outside the ring
  and the others inside it. Fixed in
  `check_interfaces` (`Piece::nonconvex_deviation`): mitred offsets of 1e-6, 1e-3, 1e-2, 0.1 and 0.25 s, each
  kept only when simple with every edge in its source edge's direction; the smallest holding the hull bounds
  the deviation by the plane distance plus its reach; the 1e-6 s offset returns at once as before, a coarser one
  is taken with the decomposition bound as min of two upper bounds. Mutant (the 1e-6 s offset only): the new
  test fails P4b. With this alone the slice-5 shaped deck (1 nm / 2 nm, s = 0.16 nm) certifies: 125304 nodes,
  753187 tetrahedra, 51001 slivers, P4 54201 faces (12493 in zones), P6 6108 pieces, 3 min 9 s Release.
- **Thin layers: pair spacing from the clearance.** The cost at small s was the pair spacing min(h, 4 s) over
  every piece, a cap the band argument does not need: a piece's pairs keep its dual faces on it while no other
  node is nearer to its plane than b = sqrt(R^2 + s^2 / 4), R <= L / sqrt 3. Each piece's clearance d (its
  distance to the pieces and boundary faces FD7's 4 s rule checks) now sets L = min(h of both regions,
  max(4 s, (sqrt 3 / 2) sqrt(d^2 - s^2))), whose band is d / 2 when above 4 s; the band is per piece (the
  generator's insertion and exemption tests and the trace pairs use their piece's). Below d = 4.73 s the spacing
  is slice 3's. All mesh fixtures still certify with their pinned deviations.
- **The decks are back on `finfet_3d.deck`'s stack (SQc as approved):** 0.5 nm SiO2, 1.5 nm HfO2, the domain
  0-85 x 0-52 x 0-14 nm, the fin z 2-12 (the shaped one tapered to z 4-10 at its top), s = 0.08 nm (FD7's 6 s on
  the ring's 0.5 nm edges). Both mesh and certify (P4, P6), Release, the two run side by side:
  - `finfet_3d_reduced.deck`: 279935 nodes, 1626836 tetrahedra, 87683 slivers, largest radius–edge 119, P4 126603
    faces (25562 in zones), P6 13146 pieces (198 in zones), 20 min 18 s.
  - `finfet_3d_shaped.deck`: 271461 nodes, 1615300 tetrahedra, 117961 slivers, largest radius–edge 2314, P4 123081
    faces (25122 in zones), P6 12832 pieces (166 in zones), 20 min 46 s.
  The shaped deck's hidden test is no longer `[!shouldfail]`. **Cost, for slice 6:** about 430 s per 1e5 nodes,
  some 14× FD15's 30 s (two meshes sharing the machine), and over twice the slice-5 meshes' nodes; FC9's "at most
  the reduced tensor mesh's nodes" and FC12 are measured there, not assumed. The thin layer's two interfaces are
  0.5 nm apart, so their spacing is about 0.43 nm (4 s would be 0.32 nm).
- Seen while measuring, not changed (both as the convex generator on `main`): `MeshReport::exempt_tetrahedra`
  counts each exemption during refinement, so it can exceed the tetrahedra; `max_radius_edge` is over every kept
  tetrahedron, not only those not exempt as its comment says (a near-degenerate one gave `inf` in a coarse
  miniature).

**Slice 6 plan, owner ruling (2026-10-30), "recommended" on SQ1–SQ5:**
- SQ1 (FC9): the generated reduced deck (about 280000 nodes at s = 0.08 nm) cannot meet "at most the reduced tensor
  mesh's nodes" (a tensor mesh whose 2× reference solves in about an hour has some 2e4 nodes). Run everything
  anyway: the reduced tensor deck, its 2× reference and the generated sweep, and record node counts and accuracy;
  FC9 stays `[!shouldfail]` under FQ2 for a ruling, its threshold unchanged.
- SQ2 (FC10): the shaped deck's full 0–0.8 V sweep at 50 mV, PARDISO on 6 threads, Release.
- SQ3 (FC8): the fixture's layers thicker than the deck's (meshable at h / 4), and s refined with h (corner zones
  are accurate to s only, so a fixed s would stall the orders).
- SQ4 (FC12, FD15): profile the generator; fix the hot spot if it is local, every mesh bit-identical; report the
  rest against FD15's budget. Larger rework is a later unit.
- SQ5 (FC8, **a change to the criterion**): as written, ohmic contacts on the fin's end faces leave the 2D
  cross-section's fin floating (refused by `Device::create`) and make the 3D field non-uniform along the axis (the
  ohmic Dirichlet value at the ends against the gated cross-section's bending), so no 2D reference exists. Instead:
  a bulk-style fin reaching the domain's bottom face, an ohmic body contact on that face (parallel to the axis, as
  the Unit 16 diode's contacts were), end faces Neumann, oxide beside the fin; the same in 2D (an edge contact).
  Gates, thresholds and the 0.5 V gate unchanged. As built the silicon is one cell, a substrate under the fin
  reaching the whole bottom face (the 2D generator refuses a patch end within 3 s of an interface junction, so a
  body patch ending at the fin's foot has no 2D reference), so the body patch has no straddled edge.

**Slice 6, second owner ruling (2026-10-30), "recommended" on SQ6–SQ7, after FC8's first measurement:**
- SQ6 (FD4, FD7, a generator defect FC8 found): `generate_complex.cpp` exempted every tetrahedron whose four
  vertices are interface nodes from refinement (inherited from the convex generator's M11, whose planar slabs never
  face each other across a thin layer). In a layer between two interfaces (HfO2, STI) the tetrahedra spanning it are
  all pair nodes and were never split, so the layer did not refine with h: an error floor (FC8's H1 and L2 orders
  about 0 on every fixture with a second insulator interface; the longest HfO2 edge 2.39 nm at h = 2 nm and 2.37 nm
  at h = 1 nm). **Fixed in slice 6**: the pair-only exemption is dropped in the complex generator only (the band
  test still protects every interface; the convex generator is unchanged, FC1); every complex mesh changes, so the
  mesh suites, the decks, FC9 and FC10 are re-run and the regression bounds re-pinned from new measurements.
- SQ7 (FC8, **a change to the criterion**): a wrapped gate bends 266° at the fin's foot, a corner where the solution
  behaves as r^0.68, so no quasi-uniform mesh reaches H1 ≥ 0.9 or L2 ≥ 1.5 (`generate2d` itself measured H1 0.86,
  L2 1.15 and 1.36 there). FC8's gate is a flat top face over an oxide fill; the interfacial oxide and HfO2 stay
  wrapped on the tapered fin, with their junctions and traces. Thresholds unchanged. The wrapped-gate case stays as
  a reported, ungated measurement.

**Slice 6, measured after SQ6 (2026-10-30; nothing committed yet):**
- FD12 pin: the convex generator's prism digest `0x5d1a9279d2a27672` (h 0.3) was read on `main` (85ad670) in a
  scratch worktree with the same check added; equal. The complex pocket digest `0x6c05e45d2f27a5df` (h 0.45) is
  pinned from this branch (FD15). The full mesh suite passes (123 cases) with the SQ6 fix and the FC12 changes below.
- FC10 (shaped deck, **before** SQ6, so to be re-run): 271,461 nodes, 1,615,300 tetrahedra, certified; equilibrium
  and 17/17 points converged; 0–0.8 V, I_D 1.33e-10 to 1.65e-5 A, swing about 62–63 mV/dec; sweep 3,866 s, total
  5,104 s. FC9's tensor parts (mesh-independent): coarse 9,367 nodes 17/17 in 25 s; reference 45,375 nodes 17/17 in
  1,007 s, I_D(0.8 V) 1.897e-5 A. FC9's generated part and FC10 are still to run on the SQ6 meshes.
- FC8, flat gate (gated), **fails** on the SQ6 meshes. Nodes 10,322 / 25,696 / 95,569 at (h, s) = (2, 0.16),
  (1, 0.08), (0.5, 0.04) nm. H1 2.64e-2, 1.96e-2, 2.02e-2 (orders 0.43, −0.04); L2 orders 1.02, 0.47; gate
  charge orders 0.91, 2.17. Wrapped gate (reported, SQ7): H1 0.55, −0.10; L2 0.87, 0.92; Q 0.85, 2.72.
  What the diagnosis has established (`generate2d` on the same fixture at (h/k, s/k), k = 1, 2, 4, 8, against the
  same reference: H1 orders 1.37, 0.72, 0.87, L2 2.28, 1.56, 1.77, so the fixture and the reference are sound):
  - The floor is a patch of about +5 mV of constant height in the 1 nm interfacial oxide around the fin's two top
    corners, its width shrinking with h; region 2 holds 15%, 28% and 76% of H1² at the three sizes. The worst 2D
    error at h = 0.5 nm is 0.69 mV.
  - It needs a permittivity contrast: with every region but the gate fill doped silicon, the orders are H1 0.66,
    0.86 and L2 1.32, 1.62 (still rising). Either contrast alone keeps it (HfO2 at eps 3.9: H1 0.44, −0.56;
    silicon at eps 3.9: H1 0.49, 0.40), both at the inner (silicon) corner.
  - Ruled out: space charge (undoped, H1 0.46, region 2 flat); dual volumes (each region's within 1e-6 of its exact
    volume, 2D and 3D); negative couplings (none); edge permittivity (the harmonic mean against the series value
    along each edge: 0.48% of the coupling weight off at h = 0.5 nm, falling with h, matching 2D near the corners);
    the slice-5 band (the slice-3 spacing min(h, 4 s) made region 2 worse, 1.29e-10 against 4.25e-11 at h = 1 nm);
    refinement points beside the corner lines (none within 2 bands of any piece edge: H1 2.02e-2 unchanged).
  - Unlike `generate2d`, which leaves the oxide sector at the corner empty (no node within 0.9 nm), the 3D
    generator puts refinement nodes there at r ≈ 0.3–0.5 nm, and the error peaks on them; removing those alone did
    not help.
- FC12 profile (temporary phase timers, removed), reduced deck at 297,758 nodes: 1,357 s, of which refinement
  1,251 s. The facet queue took 925 s: 34M pops, 1,659 splits, about 110 pops per node, almost all no-ops; then
  point location for region sizes (176–200 s, a linear scan of the complex's facets) and the band test (31 s).
  Three local changes, every measured mesh bit-identical (the mesh suite's pins; the three FC8 meshes' graph digests
  `70d9496af5218f26`, `1635249205dfbfe6`, `fa8da46bf15a8fdc` equal against a build without them):
  `Triangulation3D::tets_with_face` stops the face's tetrahedra search at the second; `ComplexIndex::locate` skips a
  facet outside its bounding box (a point on the closed facet lies in it, and a triangle's box lies in its facet's);
  `subfacet` caches each face triangle's centroid test by its local vertices in stored order (the same arithmetic).
  Refinement on the reduced deck fell to 490 s (about 600 s in all, 2.3×); FC8 at h = 0.5 nm from 136 s to about
  85 s. The rest against FD15's 30 s per 1e5 nodes (about 200 s per 1e5 on the deck): the redundant facet-queue
  entries (de-duplicating them may change the order of splits, so not provably bit-identical) and the linear
  facet scan.

**Slice 6, third owner ruling (2026-10-30), on SQ8–SQ9:**
- SQ8 (FC8): keep investigating the corner error before the rest of the slice; a generator fix comes back for a
  ruling.
- SQ9 (FC12): keep the three bit-identical changes, confirm the digests on the mesh suites and decks, and report
  the remainder against FD15 for a later unit.

**Review of slices 1–6 (2026-10-30, owner "fix"; nothing committed yet):** a code review of the branch against
`main` found nine items; each is fixed or measured as below. Every change leaves the meshes bit-identical: the mesh
suite's pins, seven interface fixtures and an 88,835-node gated liner fin (graph digests equal against a build
without the change), and the three FC8 meshes (`70d9496af5218f26`, `1635249205dfbfe6`, `fa8da46bf15a8fdc`, as
recorded above). Release 127/127.
- **Deck stations matched exactly** (`build.cpp`): a range end spelled `300nm` against a station `0.3um` was "not a
  station". The range ends now take the nearest station within 1e-12 of the stations' scale; `mesh::extrude`'s
  patch boxes hold vertices within 1e-12 of the extrusion's extent. New deck test (lowest spelling for the station,
  highest for a box face); both exact-comparison mutants fail it.
- **Extraction's centroid test** (`generate_complex.cpp`): it summed `q/3` per vertex while refinement's
  `subfacet` used `(a+b+c)/3`, so a centroid within rounding of a face's boundary could be judged differently.
  Extraction now asks `subfacet` for any triangle of the face triangulation (the old test remains for others). No
  test: no fixture puts a centroid within an ulp of a face edge.
- **The centroid cache grew without bound** (an `unordered_map` keyed by packed vertex ids, uncached above 2^21
  vertices): now a memo per triangle slot holding the vertices it was taken with.
- **The run record of a reloaded mesh** said "file" next to an extrude3d generation digest: a file mesh with a
  digest now records `extrude3d`, the generator the digest belongs to (FD14 deck test extended).
- **GUI vocabulary**: `patch` in an extrude3d `[mesh]` offered generate2d's `<x y>, <x y>`; it now has its own
  entry `<x0 y0 z0 x1 y1 z1>` (new test); the `radius_edge` menu text names extrude3d.
- **Duplicated code**: one walk around a vertex (`Triangulation3D::walk_around`) now serves `incident`,
  `tets_with_face` and the new allocation-free `tets_with_edge` (`segment_encroached`). The generator uses
  `polygon2d.hpp`'s exact `collinear` and `in_triangle_2d`. Its `circumcentre_fast` and `min_dihedral` are **not**
  copies of `generate3d`'s (other arithmetic order and formula); sharing them would change every mesh, so they stay.
- **The volume-overlap check** copied every vertex of the complex per pair of volumes: `ComplexIndex::locate` now
  takes a facet list (private), the same tests in the same order.
- **Facet-queue de-duplication** (a "queued" flag per face triangle slot): built and measured, bit-identical on
  every mesh above, but no faster (88,835 nodes: 26.6 / 24.0 s with it, 24.1 / 22.7 s without), so not kept.
  The redundant pops are cheap since the FC12 caches; the deck's remaining time is elsewhere (not yet profiled).

**Owner ruling on FC8 (2026-10-31), "skip fc8":** the investigation (SQ8) stops. FC8's flat-gate case fails and
is an **accepted limitation**: the test stays `[!shouldfail]` with its thresholds unchanged and its measurement
above (H1 orders 0.43, −0.04; L2 1.02, 0.47; gate charge 0.91, 2.17 at (h, s) = (2, 0.16), (1, 0.08),
(0.5, 0.04) nm). What is known: a constant-height error of about 5 mV in the 1 nm interfacial oxide at the fin's
two top corners, present only with a permittivity contrast; `generate2d` on the same fixture converges. The
wrapped-gate case stays reported (SQ7). A later unit may take it up; the cause is not found.

**Slice 6, completed (2026-10-31):**
- **FC13, hand mutants** (Release, a separate build tree; each mutant built and run against its suite's tests):
  - Before this session (2026-10-31): FD1, FD2, FD3, FD6 and FD7's sector rule killed; FD5 (the node-cap test
    stopped by the corner count, `max_nodes` 20 below the L-prism's 24 vertices), FD7 (4 s loosened to 3 s: the
    thin-layer and near-boundary cases were refused first by the crossing rule), FD11 (a region's h) and FD12/FD15
    (no pins) survived; each closed by a new or corrected test (`b4667ff`; FD12 and FD15 digest pins `019da1f`).
  - FD13 (deck): the duplicate (cell, section) check, t_lo < t_hi, a patch box's lo <= hi, the station tolerance,
    "needs at least one cell": five mutants. Three survived (`t_lo == t_hi`, the box's z only, no cells); new
    deck tests close them; all five killed.
  - FD14 (digest, io, run record): the volume region and radius–edge dropped from the generation digest, the
    mesh file's digest record, the deck's `extrude3d` type name for a file mesh, the run record's digest: all five
    killed by the existing tests.
  - FD4 (extraction): the patch-of-node fallback (a node in two patches' triangles) killed. **Two survivors,
    argued**: the flood fill's skip of flat tetrahedra and `is_subfacet`'s rejection of a triangle with three
    vertices on one complex edge. Both concern zero-volume tetrahedra (all four vertices on one face, or three on
    one edge). A flat tetrahedron's triangles lie in its face's plane; where the face is convex around them, each
    is a subfacet (its centroid is in the face), so the flood cannot cross it with or without the guard. Only a
    flat tetrahedron spanning a reflex notch of a non-convex face could let the flood through, and then interior
    tetrahedra are lost and the volume check (kept volume against the complex's within 1e-12) refuses the mesh
    (`degenerate_mesh`): without the guards a difference is a refusal, never a wrong mesh. No fixture reaches
    them; kept as defensive guards.
  - FD8, FD9, FD10: mutated in slices 3–4 (above), each killed.
- **FC12, cost** (Release; new `[.benchmark]` test on the gated liner fin, measured while FC10 ran on 6 threads):
  14,562 nodes in 4.00 s and 88,835 nodes in 23.89 s, about 27 s per 1e5 nodes, **within FD15's 30 s**; exact
  location of every node against the whole complex (FD6) 0.15 s and 1.11 s. The FinFET decks are slower per node
  (about 130 and 160 s per 1e5: the reduced deck built in 384 s at 297,758 nodes, the shaped in 454 s at 288,729, mesh, certificate and device together): their thin layers at s = 0.08 nm put most nodes in interface bands, where the facet
  queue and the per-node region size dominate (the SQ9 profile above); a later unit.
- **SQ9 on the decks** (the three FC12 speedups and the review fixes bit-identical): a scratch build of `71e6e98`
  with only SQ6 applied (and a digest print) against this branch's head:
  the reduced deck `d5a99737270489bf` (297,758 nodes) and the shaped deck `719a79d9189eba8a` (288,729 nodes) on
  both: **bit-identical**. Pinned in a new hidden test (`[.fc_digests]`, `tests/deck/complex_finfet_test.cpp`). The
  pre-speedup build meshed each deck in about 21 min (sharing the machine with FC10).
- **FC10 on the SQ6 meshes** (`[.fc10]`, Release, PARDISO on 6 threads): **met**. 288,729 nodes, 1,695,931 tetrahedra, 45,214
  slivers, certified (P4 123,121 cross faces, P6 12,968 patch-boundary pieces); 17 of 17 points converged, one bias
  sub-step; build 454 s, sweep 3,099 s. I_D 1.94e-10 A at 0 V to 1.37e-5 A at 0.8 V; about 62 mV/dec from 0 to
  0.15 V. Against the pre-SQ6 run (1.33e-10 to 1.65e-5 A) the on-current is 17% lower and the off-current 46% higher:
  the thin layers now refine (SQ6), which changes the mesh everywhere near the gate stack. No reference exists for
  this shape (no accuracy gate).
- **FC9 on the SQ6 meshes** (`[.fc9][!shouldfail]`, Release, PARDISO on 6 threads): **fails, as pre-stated (FQ2)**. Every
  run converged (17/17 each). Coarse tensor 9,367 nodes (sweep 14 s), tensor reference 45,375 nodes (488 s),
  generated 297,758 nodes (build 384 s, sweep 2,190 s). Generated against the reference: I_D at 0.8 V 1.634e-5
  against 1.897e-5 A (−13.9%; the coarse tensor +1.2%), at 0 V 2.56e-10 against 2.18e-10 A (+17%); V_T 0.200 against
  0.231 V; max g_m 5.48e-5 against 7.63e-5 A/V (−28%); the swing 64.4 against 64.9 mV/dec, the only quantity within
  its bound. The node rule fails by 32×. The generated mesh is further from the reference than the tensor mesh
  with 1/32 of its nodes: a systematic shift (lower drive, earlier turn-on), not a resolution floor. Candidates,
  not diagnosed: the missing structured layer under the gate (item 5, the risk stated before the run) and the
  FC8 corner error at the fin's top corners. **For an owner ruling** (FQ2: pull item 5 in, or change the gate; the
  threshold is not loosened).
- **FC14**: **not met; open for an owner ruling.** Re-measured serially in Debug (2026-10-31): `[complex_generate]`
  38.0 s, `[complex_interfaces]` 67.0 s, `[complex_patches]` 17.6 s, `[complex]` 7.9 s, the deck's `[extrude]`
  6.4 s, and `mesh.interfaces` about 9 s more than before the unit: about 146 s against the 120 s budget. Slice 4's
  "about 102 s" counted slices 3 and 4 only; slices 1 and 2 (`mesh.complex`, `mesh.generate_complex`, about 50 s)
  were never added in. Cut in this slice without losing a check: the FC13 closed-shell fixture built from nine cells
  instead of a 5 × 5 grid (13.2 s to 3.6 s Debug; the 4 s → 3 s mutant still killed), and the two determinism
  cases mesh once against a pinned digest instead of twice (the pocket `0x6c05e45d2f27a5df`, the liner fin
  `0xfffed114cf3222c1`, new). The remaining 26 s needs a gate case moved to `[.benchmark]` (then run in neither
  ctest suite): the largest are the gated section with junction points and non-convex pieces (FC6, FD9; 22.7 s),
  the liner fin (FC6; about 7 s now), FC5's linear field (13.3 s), the tapered ring (8.6 s) and the trace-crossing
  patch corner (9.6 s). Or the budget is raised.
- **Suites** (2026-10-31, `/W4 /WX`, 6 jobs, before the FC14 cuts): Release 127/127 (804 s), Debug 116/116
  (1,071 s). After the cuts, the changed mesh tests pass in Release (18 cases) and Debug.

**Owner ruling on FC9 and FC14 (2026-10-31), "skip FC9 and FC14":** both accepted as limitations. FC9 stays
`[!shouldfail]` with its measurement and threshold unchanged; FC14's Debug growth stays at about 146 s, no case moved
to `[.benchmark]`. Unit 53 is ready for review with FC8, FC9 and FC14 unmet and accepted.
