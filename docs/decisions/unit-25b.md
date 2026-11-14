<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 25b — derived fields, isosurfaces and streamlines (proposed 2026-10-21; on `render/derived-fields`)

**Owner request:** VQ7 (2026-10-17): "25b, derived |E| and |J|, isosurfaces and streamlines", after Unit 45
(VQ12); gates VG4 (in 25b also |E|, |J| with streamlines and isosurfaces) and VG8 (scientific validation), with VG6
(suites and mutants). The owner's "GO" (2026-10-21) asked for this proposal. Documentation only until the owner
approves.

**What exists (audited on `main` at `81622ca`).**
- Every result with fields holds what the vectors need, so nothing is solved again: the potential per node, and the
  electron and hole current through every mesh edge (`BiasPoint::edge_current_n/p` and
  `TransientSnapshot::edge_current_n/p`, in every result file since Unit 9). The mesh gives each edge its length and
  its coupling area (`mesh::Edge::coupling_area`), so I_ij / A_ij is the current density along the edge, in A/cm² in
  every dimension.
- `render::Geometry` has the linear elements (Kuhn tetrahedra of tensor cells, generated simplices), the slices, the
  barycentric interpolants and marching triangles (2D contours). It has no marching tetrahedra, no point location
  and no vector field.
- `render::Field` is nodal; `field_values` returns one value per node; `nitcad plot` draws field, slice and contour
  views of it. Line batches have one color and a pixel width.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| DS1 | **E per element:** E = −∇ψ of the linear interpolant on each segment, triangle or tetrahedron, in V/cm; insulators included. | Exact for a linear potential (VG8); the field in the oxide matters as much as in the silicon. | Constant per element; an element whose nodes lie in two regions (the box method puts the interface between nodes) mixes both permittivities' fields. |
| DS2 | **J per element from the stored edge currents:** for each element, the vector whose projections on its mesh edges best match the edge densities I_ij / A_ij (least squares over the edges with a positive coupling area and both ends in a semiconductor); Jn and Jp separately, J = Jn + Jp. An element without D independent such edges has no current vector. | The solver's own currents (Scharfetter–Gummel), not a second discretization; exact for a uniform current (VG8). A Kuhn tetrahedron of a tensor cell has exactly three axis edges, so the fit is a solve there. | Conduction current only: a transient snapshot's displacement current is not stored (named in the title). |
| DS3 | **Nodal values for drawing:** each node gets the average of its elements' vectors, weighted by the volume of its control-volume piece in each element (the region rule of 25a); \|E\|, \|J\|, \|Jn\|, \|Jp\| are the magnitudes of those averages. New `render::Field` keys `E` (V/cm), `J`, `Jn`, `Jp` (A/cm²), log scale by default; 1D gives profiles. | Every 25a view (maps, 3D surfaces, slices, contours) draws them unchanged; the average stays inside the node's own region. | Smoothed: a peak narrower than an element is spread (the element values remain available to the code and the tests). |
| DS4 | **Isosurfaces (3D):** marching tetrahedra on the elements, for any nodal field, one or more levels; triangles oriented toward increasing value, shared edges cut once (watertight), each level colored by its colormap position and drawn translucent over the solid. Values are interpolated in the field's scale (log10 for densities and \|J\|), as the colormaps already are. In 2D the same request draws contour lines (25a). | Isosurfaces are the standard 3D view of depletion edges, channels and current paths; interpolating in the drawn scale puts a 1e16 surface where the map shows 1e16. | Exact planes on a linear field only on a linear scale. |
| DS5 | **Streamlines of E, J, Jn or Jp (2D and 3D):** the nodal vectors interpolated linearly in each element; adaptive Runge–Kutta (Dormand–Prince 5(4)) with point location by walking element neighbours from a bucket grid; a line stops at the domain boundary, a contact, the edge of the current-carrying region (J), a vanishing field or a length limit. Drawn as lines, in one color over the field map, or over the translucent solid in 3D. | Integrating the same linear field the maps show keeps lines and colors consistent; an adaptive integrator gives VG8's closure test a stated tolerance. | Lines have no arrowheads (vector glyphs are 25c). |
| DS6 | **Seeds:** J lines start on a contact, spaced so each carries an equal share of its current (the line density then shows the current density); E lines start evenly along a contact or gate; an explicit seed line or plane (`--seeds`) for either. | What a TCAD user reads from a current-line plot is where the current flows and how much. | Equal-current seeding needs the contact's current density along it, read from its edges. |
| DS7 | **`nitcad plot`:** `--field E|J|Jn|Jp`; `--streamlines E|J|Jn|Jp` over any field view (`--lines <n>`, `--from <contact>`, `--seeds <line or plane>`); `--view iso --field <f> --levels <v,...>`. | One more switch per feature on the existing command. | — |
| DS8 | **Out of scope:** vector glyphs, clipping, exploded views and filmstrips (25c); displacement current; per-element drawing (DS3's alternative); animated or tube-shaped streamlines; E or J stored in results (the `render` layer reconstructs them). | Keeps the unit to VQ7's list. | — |

**Owner decisions (asked):**

| # | Question | Recommendation |
|---|---|---|
| BQ1 | Draw \|E\| and \|J\| as nodal averages (smooth, DS3) or as one value per element (faceted, exact)? | **Nodal averages**, as commercial viewers draw them; the exact element values are what the gates test. |
| BQ2 | Seed J streamlines on a contact with an equal current per line (DS6), or only evenly along a user's line? | **Equal current per line from a contact**, with `--seeds` for an explicit line: the line density then means current density. |
| BQ3 | Interpolate isosurfaces in the field's scale (log10 for densities and \|J\|), or always linearly in the value? | **In the field's scale**, so the surface lies where the map's color says the level is. |

**Acceptance gates, fixed before any implementation (VG8, VG4 and VG6 made precise):**

| # | Gate |
|---|---|
| BG1 | **E and J exact per element (VG8):** \|E\| of a linear potential is exact (to rounding) on every element of 1D, 2D and 3D tensor and generated meshes; J from edge currents of a uniform current (I_ij = J · t_ij A_ij) is exact on every element; a solved uniform resistor bar gives \|J\| = I / area at every node to the solver's current resolution; Jn + Jp = J. |
| BG2 | **Isosurfaces (VG8):** marching tetrahedra on a linear field give the exact plane (every vertex on it, normals along the gradient); on a sphere's distance field inside the mesh the surface is closed (no edge used once) and outward, and its enclosed volume converges at second order under refinement (rate measured on three meshes, between 1.8 and 2.2). |
| BG3 | **Streamlines (VG8):** a uniform field gives straight lines to rounding; the rotational field (−y, x) (linear, so exactly represented) closes after one turn within the integrator's tolerance; lines stop at the boundary and at contacts; equal-current seeds of a uniform bar are evenly spaced, and of a known nonuniform contact current spaced by its integral. |
| BG4 | **Images (VG4):** WARP references, kept like 25a's: the 2D MOSFET's \|E\| with E lines and its \|J\| with current lines; the 3D MOSFET's isosurfaces of n and its current lines over the translucent solid; within 25a's per-pixel tolerance; the D3D12 debug layer clean in Debug. |
| BG5 | **No numerical change:** nothing outside `render` and `nitcad plot` changes (a header-boundary test as VG1); the assembly traces byte-identical to `main`'s; every suite passes unchanged. |
| BG6 | **`nitcad plot`:** each new field, the isosurface view and streamlines written from decks and from result files; a result without edge currents, a 2D isosurface request and an unknown contact refused with a message (exit 1). |
| BG7 | **Benchmarks:** \|E\| and \|J\| reconstruction, isosurfaces and 100 streamlines on the 3D MOSFET and on a 10⁶-tetrahedron mesh, Release, recorded. |
| BG8 | **Suites and mutants (VG6):** Release and Debug at every commit, no warnings under `/W4 /WX`; mutants in the gradient, the edge fit, the averaging weights, marching tetrahedra (an orientation, a shared-edge cut), the integrator (a stage weight, the step control), point location and the seeding, all caught. |
| BG9 | **No scope creep:** nothing outside DS1–DS7; ARCHITECTURE.md 6.14 and 11.3 as built after acceptance. |

**Commits (proposed):** 1) element and nodal E and J, the four field keys and `--field` (BG1, BG5, part of BG6);
2) marching tetrahedra and `--view iso` (BG2); 3) point location, the integrator, seeding and `--streamlines`
(BG3); 4) the references and the benchmark (BG4, BG7). About 1,700 lines with tests (VQ7's estimate).

**Risks:** element-mixed interfaces blur E across an oxide–silicon interface on coarse meshes (DS1; the gates use
single-region elements, the references show it as drawn); streamlines through regions of near-zero current wander
(stopped by DS5's vanishing-field rule); point location on large 3D meshes (BG7 measures it).

**Status:** authorized by the owner 2026-10-22 ("Approve, as recommended"): BQ1 \|E\| and \|J\| drawn as nodal
averages, BQ2 J streamlines seeded on a contact with an equal current per line (and `--seeds` for an explicit line),
BQ3 isosurfaces interpolated in the field's scale.

**Commit 1 (\|E\| and \|J\|; BG1, BG5, BG6 in part), as built.** `vectors.hpp`: `element_vectors` (E = −∇ψ per
element, by Cramer's rule on its edge vectors; Jn and Jp by least squares on the element's mesh edges with a positive
coupling area and both ends in a semiconductor, J their sum), `control_volume_shares` (tensor: binomial(D, m) / 2^D
for the vertex m steps along the Kuhn path, which is the corner sub-box's part of the simplex, so the shares add up to
the box method's control volumes exactly; simplex: 1 / (D + 1)), `nodal_vectors`, `magnitudes`. `render::Field` gains
`E`, `J`, `Jn`, `Jp` (log10), given by a `field_values` overload taking the geometry (the old one refuses them).
`nitcad plot --field E|J|Jn|Jp`: on a log scale without `--range`, the ten decades below the largest value, values
under them at the bottom color (below that the reconstructed currents are rounding: Jn at equilibrium spans 1e-42 to
1e-10 A/cm²); a magnitude that is zero everywhere on a linear 0..1 scale.
- **BG1:** E of a linear potential and J from the edge currents of a uniform current exact (10⁻¹¹) on every element
  and node of 1D, 2D and 3D tensor and 2D and 3D generated meshes; the tensor shares sum to the mesh's control
  volumes (10⁻¹²) in 1D, 2D and 3D; a nonuniform 1D E averages each node's two half segments; a solved 2D resistor
  bar at 0.1 V has \|J\| = I / height at every node within 10⁻⁶ plus the current resolution, along +x; Jn + Jp = J;
  elements with an oxide node have no current, oxide nodes NaN \|J\|; results without edge currents and other or
  element-less geometries refused.
- **Finding (DS2 amended):** a quasi-static result stores zero edge currents, not none (`sweep_bias` reports zero
  currents), so its \|J\| draws as zero; the refusal is for results without edge currents.
- **BG6 (in part):** E, J, Jn and Jp drawn from a 2D deck and from its result file (the equilibrium point's \|J\|,
  exactly zero, on 0..1), a 1D \|E\| profile (the junction's peak), a quasi-static \|J\| drawn as zero.
- **Mutation checks (BG8, Release): 14 of 14 caught** — E's sign, the 2D and 3D Cramer solves, insulator edges kept
  at either end, a reversed edge tangent, no division by the coupling area, holes subtracted, equal tensor shares, a
  wrong Kuhn step count, unweighted volumes, a 2D magnitude, `nitcad plot` on the old overload, the zero-field
  scale off. (With oxide only at the top, the first-node check was unseen: tensor numbering puts the silicon node
  first on every such edge, so oxide was added below. Two first versions failed only the `/WX` build, an unused
  variable, and were replaced by mutants the tests catch.)

**Commit 2 (isosurfaces; BG2, BG6 in part), as built.** `geometry.hpp`: `isosurface` (marching tetrahedra into an
indexed mesh: one vertex per cut mesh edge, computed from its lower-numbered node so every element around the edge
makes the same point, or the node itself for a crossing exactly at it; a node at the level counts as above it, as in
the 2D contours; triangles of zero area dropped; each triangle turned toward higher values by the sign of
Σ (v_k − level)(x_k − p₀), which for a linear field is the gradient's). `ColorScale::transform` is public.
`scene.hpp`: `iso_view` (3D: the surfaces of the field's scale coordinate (BQ3), colored by each level's colormap
position, smooth-shaded, opaque for one level and translucent for more, in the device's faint outline with the
contacts; 2D: the field map with contour lines at the levels). `nitcad plot --view iso --field <f> --levels
<v>[,<v>...]`.
- **BG2:** a linear field's surfaces on graded tensor and generated 3D meshes: every vertex on the plane (10⁻¹²),
  every normal along the gradient (cosine ≥ 1 − 10⁻⁹), no directed edge twice. A sphere's distance field (R = 0.3 µm,
  centre off the grid planes) on 8³, 16³ and 32³ cells: closed and consistently oriented, outward, volume errors
  8.81e-2, 2.15e-2 and 5.41e-3, **rates 2.03 and 1.99**. Crossings at nodes give one vertex per node (both edge
  orientations); a level at the field's maximum gives the boundary face; NaN nodes and 2D geometries give nothing.
- **BQ3:** a field exponential in y (log scale) has its 10^14.3 surface on the plane y = 0.43 µm between grid
  planes, and its 2D contour on that line.
- **BG6 (in part):** 3D surfaces of n (two levels) and ψ (one) from a deck, 2D contours at a level; refused: no
  `--levels` (usage), an unreadable level (usage), a negative level on a log scale, a 1D device. **For the owner:**
  BG6 listed "a 2D isosurface request" among the refusals while DS4 draws contours at the levels in 2D; DS4 is built
  (2D draws contours), and the refusal applies to 1D.
- **Mutation checks (BG8, Release): 9 of 9 caught** — the shared edge cut from either end, the node rule off, a node
  at the level counted below, a crossed quadrilateral, flipped orientation, contour levels ignored, interpolation in
  the value instead of the scale, one level translucent, `nitcad plot` allowing 1D. (The first two survived a first
  run: Kuhn simplices number their nodes in increasing x, so only a decreasing field reaches the node rule's other
  branch, and only a level at the field's maximum shows the at-level convention; both cases were added.)

**Commit 3 (streamlines; BG3, BG6 in part), as built.** `streamlines.hpp`: `Locator` (the hint element, then the
elements sharing a node with it, then a bucket grid of about one element per bucket over the elements' bounding
boxes; barycentric weights by Cramer's rule, inside at −10⁻¹⁰); `trace_streamline` (the nodal vectors interpolated
linearly, followed at unit speed, adaptive Dormand–Prince 5(4) with first-same-as-last, local error per step at most
`tolerance` (10⁻⁵) times the device's size, steps at most half the local element's size so the drawn polyline follows
the field; forward and backward from the seed; stops at the boundary (to 10⁻⁷ of the size), at an undefined or
vanishing field (|v| ≤ 10⁻¹² of the largest), at a sink (a step that turns the direction back), at 20 device sizes
or 20,000 steps; |v| kept per point); `contact_node_currents` (the stored edge currents from each contact node to
nodes outside the contact); `contact_seeds` (BQ2: the contact's boundary segments (2D, walked in order along each
chain) or faces (3D), each node's \|current\| spread evenly over its part, seeds at the midpoints of equal shares;
without currents, by length or area); `segment_seeds`, `plane_seeds` (by area over `Geometry::slice`).
`scene.hpp`: `add_streamlines` (one color, or each segment through the colorbar in 64 bins), `outline_view` (the
outline shell, now shared with the slice and isosurface views, and the contacts). `nitcad plot --streamlines
E|J|Jn|Jp [--lines n] [--from contact] [--seeds a:b | x=um]`: 2D over the field map (by default the lines' own
magnitude) in near-black; 3D over its own magnitude, the lines colored by it in the device's outline, else over the
chosen view. Default seeds: J lines on the contact carrying the most current, E lines on the first gate or electrode,
else the most biased contact.
- **BG3:** a uniform field's lines straight to rounding, ordered, ending on y = 0 and y = 1 µm (10⁻⁶ of the size);
  the rotational field (−y, x) closes after one turn to 1.2e-15 cm (109 points) and, on a coarse mesh at tolerance
  10⁻¹², to 1.4e-16 cm; a field with a kink on a node column, v = (1, k\|x − xm\|), follows its closed-form line within
  4.2e-14 cm (budget: tolerance × size × points); lines stop at a NaN region (no field), at a sink, and a zero field
  gives the seed alone; 2000 points located in a generated 3D mesh, their weights rebuilding them to 10⁻¹²; equal-
  current seeds of a uniform current at (i + ½)/N of the contact, and of j(y) = 1 + 3y/H at equal shares of its
  integral (10⁻⁹), denser where the current is; the two contacts' currents balance (10⁻¹²); plane seeds on the plane.
- **BG6 (in part):** J, Jn (from a named contact) and E (from a segment, over ψ) lines in 2D; J lines, E lines from a
  plane and Jp lines over a slice in 3D; refused: an unknown field, a malformed segment, a 3D segment in 2D, a plane
  in 2D, `--lines 0` (usage), an unknown contact, a plane that misses the device, a 1D device.
- **Mutation checks (BG8, Release): 12 of 12 caught** — the quantiles' midpoints, a Dormand–Prince stage weight and
  the error estimate's last weight, the barycentric weight of node 0, the length limit, every step accepted, the
  sink stop, the step cap, the backward part's order, a contact current's sign, the density not per length, a half
  segment reversed. (Every step accepted survived until the kinked field: in a smooth field the step sizing keeps
  each attempt within the tolerance, and only a kink makes a step fail.)

**Commit 4 (references and the benchmark; BG4, BG7), as built.** The four references showed two faults of the lines
the unit tests had not, both fixed here in `vectors.hpp`'s `streamline_field`, which the plot and the references now
trace instead of the bare nodal vectors:
- **Finding 1, the elements along an interface:** with insulator nodes NaN, every element joining the last silicon
  row to the first oxide row has no current vector, so each current line stopped at the silicon surface, short of the
  drain contact (in an on-state MOSFET it would cut every channel line). Insulator nodes now carry zero current for
  streamlines (the |J| map keeps them as no data), so lines continue along the interface and stop only where the
  current vanishes in the oxide.
- **Finding 2, normal components on no-flux boundaries:** the element average at a boundary node comes from one side,
  so it kept a normal part there (the top silicon row's current pointed into the oxide, and lines left through it
  0.25 µm short of the drain). The model has no normal current at semiconductor–insulator interfaces and no normal
  current or field on the outer boundary away from contacts (homogeneous Neumann), so `streamline_field` removes the
  normal part there: per node its boundary pieces' normals, grouped within 30° (one face, or a gently curved one)
  and orthonormalized (the faces meeting at an edge or corner each count), every component along them removed;
  contact nodes keep theirs. DS5 is amended accordingly.
- **BG4:** WARP references, 960 × 720, in 25a's tolerance (at most 0.5% of pixels off by more than 24):
  `mosfet2d_efield_lines` (\|E\|, 16 lines from the gate through the oxide to the depletion edge, where they end at
  the field's sink), `mosfet2d_current_lines` (\|J\|, 16 equal-current lines from the source: this subthreshold
  device conducts under the channel through the 1e15 substrate, then along the drain pocket's surface to the drain),
  `mosfet3d_iso_n` (n = 1e14 and 1e18, the source and drain junctions), `mosfet3d_current_lines` (24 lines colored
  by \|J\|, source to drain contact). Each was looked at before being kept; all 15 references pass; no D3D12 debug
  errors.
- **Tests:** insulator nodes zero and interface nodes tangential with J = (5, 2), and a line along the top silicon
  row reaching the far contact inside the interface band; E = (3, 4) on 2D and 3D bars keeps (3, 0) on the
  non-contact sides, edges and corners included, and (3, 4) at contacts and inside.
- **BG7 (Release, other jobs running, so upper bounds):** the 3D MOSFET (3,087 nodes, 14,400 tetrahedra): E 0.001 s,
  J 0.003 s, an isosurface 0.0004 s, streamline field and locator 0.009 s, 100 lines 0.012 s; a 3D tensor device of
  185,193 nodes and **1,053,696 tetrahedra**: E 0.048 s, J 0.495 s, an isosurface 0.014 s (41,078 triangles),
  streamline field and locator 0.73 s, 100 lines 0.054 s (24,844 points).
- **Mutation checks (BG8, Release): 5 of 5 caught** — insulator nodes not zeroed, every normal merged into one,
  interfaces ignored, no projection, contacts projected. Unit 25b in all: 40 of 40.

**BG5 (no numerical change), as checked.** Only `render`, `nitcad plot`, tests and documents differ from `main`
(`81622ca`); `CMakeLists.txt` adds the three sources and test files. Assembly traces of the bias, sweep and transient
tests (`--order decl`): the 38 cases other than "bias: assembly cost on the legacy diode" are byte-identical
(666,691 lines); that case repeats residual evaluations for a time budget, and on the branch (run while three test
suites shared the CPU) it stopped after 118,371 of main's 185,955 lines, which are an exact prefix of main's. Its count
depends on the machine's load, not on the code.

**BG8 suites:** Release 22/22 and Debug 16/16 at each of the four commits (commits 2–4 in clean worktrees), no
warnings under `/W4 /WX`; the Debug render suite includes the D3D12 debug layer check (no errors).

**Acceptance (2026-10-22).** The owner pre-accepted the unit and authorized push, PR and merge, then the start of 25c
without a further approval ("As Unit 25b is done, accept it push, open pr and merge, then start 25c without my
approval"), given while only commit 1 existed; it was taken as acceptance on the gates passing, which they did.
ARCHITECTURE.md 6.14, 11 and 11.3 as built (BG9). For the owner's record: DS2 amended (a quasi-static result's
currents are zero, not absent), DS5 amended (the streamline field, commit 4), and BG6's "2D isosurface request
refused" superseded by DS4 (2D draws contours at the levels).
