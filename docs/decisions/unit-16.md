<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 16 — unstructured 2D and 3D meshing (branch `mesh/unstructured`, revised 2026-10-14; on `main` as `a09cf43`, PR #35)

**Owner rulings (2026-10-14), unchanged:** R1 2D and 3D in this unit; R2 an in-tree mesh generator, no new dependency;
R3 material interfaces as node regions with straddling node pairs; R4 analytic doping profiles included. The branch
starts at `main` `23199a9`. Everything below is proposed and **awaits the owner's ruling before any code**. Revision 2
follows the review of revision 1 (`cae0ef0`): its M2 contradicted M3, its sliver refusal contradicted M3, it had no
tie-break rule, and M7 lacked the conditions under which it is exact.

**Evidence labels.** *Proven*: follows from a stated argument or a cited theorem whose hypotheses are checked on every
generated mesh. *Cited, to verify*: from the literature; the theorem's exact hypotheses must be confirmed against the paper
before `ARCHITECTURE.md` cites it. *Heuristic*: no guarantee; bounded by a cap and checked after the fact.
*Limitation*: known and stated, not fixed here.

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| M1 | Box method on the circumcentric dual. An edge's coupling area is the signed sum, over the simplices sharing it, of the dual-face pieces. A node's volume and a boundary node's face area are the signed sums of its pieces. Circumcentres are computed with exact determinants (Shewchuk expansion arithmetic) and one final rounding. | Proven: on a Delaunay mesh the dual is the Voronoi diagram, which depends on the points only, not on which tetrahedra were chosen. | Standard Voronoi finite volume; reduces to the tensor-grid box method on split rectangles and Kuhn-split boxes. Exact determinants keep each circumcentre accurate to a few ulps however flat the element. | — |
| M2 | Checks on **totals**, never per-element pieces (an obtuse triangle's pieces are legitimately negative). Tolerance τ = 1e-12, scale L_e^(D−1) for an edge and the node's (local mean edge)^D for a volume. A total within τ·scale of zero is zero: the edge is dropped, a zero volume is refused. A total below −τ·scale is refused (`degenerate_mesh`, index of the edge or node). Roundoff of order 1e-17 is never refused. | Proven that totals are ≥ 0 under M3's conditions; τ is about 4.5e3 ulps, the precedent of `length_tolerance` in `mesh.cpp`. | Fixes revision 1's contradiction with M3; cocircular diagonals (exact zeros) are dropped, not refused. | A total that is truly negative but within τ of zero is treated as zero. |
| M3 | Guarantees. (a) Delaunay ⇒ every interior total coupling ≥ 0. (b) Delaunay plus the Gabriel property for every boundary segment or facet ("boundary conforming Delaunay") ⇒ the dual is the Voronoi diagram clipped to the domain: every volume > 0 and every boundary area ≥ 0. (c) The generator checks (a) and (b) exactly on its output (P1, P2 below), so they are postconditions, not just claims. | (a) proven; (b) cited, to verify (Gärtner, Si, Fuhrmann 2010; Fuhrmann's lecture notes). | (b) is the property Voronoi box methods require. | — |
| M4 | Exact predicates (orient2d/3d, incircle, insphere) ported from Shewchuk's public-domain adaptive code. Ties (incircle or insphere = 0) are broken by the symbolic perturbation of Devillers and Teillaud (2011): perturb the lifting, decide by orientation tests on subsets ordered by node index. orient tests are never perturbed, so no tetrahedron of four coplanar points can exist; every element has orient > 0 (P0). | Proven for the predicates and for determinism; the tie-break is cited, to verify. | Tensor lattices, Kuhn splits and mirrored pairs are exactly cocircular or cospherical by design; ties must resolve the same way on every machine (6.8). | Each tie costs a few extra orientation tests. |
| M5 | 2D domains: any planar straight-line graph (non-convex polygons, holes) by constrained Delaunay triangulation, then refinement until no segment is encroached (so the CDT is Delaunay). | Cited, to verify: boundary conforming Delaunay meshes exist and are computable for arbitrary polygonal regions in 2D (Gärtner–Si–Fuhrmann abstract). | Real cross-sections. | — |
| M6 | 3D domains: convex polyhedra only. Boundary conformity (M3 b) is guaranteed only when every input angle is at least arccos(1/3) ≈ 70.53°; boxes (90°) satisfy this. Convex domains with smaller dihedral angles are allowed and handled by the heuristic in M11, with P2 still checked. Non-convex domains are refused (`invalid_input`) and deferred. | Proven: for convex domains the Delaunay tetrahedralization meets the hull without boundary recovery. 70.53° is cited and checked (Gärtner–Si–Fuhrmann 2010, U6). | No 3D constrained boundary recovery (the TetGen-class problem). Interfaces need no conforming faces (M7). | The angle is the theorem's hypothesis, not a refusal rule. |
| M7 | (Ruled R3.) Every node has one region; nodes are never on an interface. On each interface, a pair set: a 1D (2D domain) or 2D (3D domain) point set P on the interface with covering radius R(x), and for each p ∈ P the nodes p ± (s/2)·n. A **protected band**: no node outside the pair set lies within b(x) = √(R² + s²/4) of the interface. Pair spacing s = min(h(x), lfs(x)/2), where lfs is the local feature size, so pairs from different features cannot interfere. | Proven (away from corner zones, M8): with mirror symmetry and the band, every point of the interface is nearest to a pair node and every off-interface point is strictly nearer its own side's mirror image. So the Voronoi faces between the two regions lie in the interface, couplings between non-mirror nodes across it are exactly 0, and each cross edge's midpoint lies on the interface (the assembler's assumption, `assemble/interface_edges.hpp`). | Generalizes `mesh::straddle_interface` exactly. | Coordinates p ± s/2 round, so cross couplings come out at about 1e-17 relative instead of 0 (absorbed by M2's τ). The interface is a polyline or polygon at resolution s. |
| M8 | Corners and junctions: a bend in an interface; an interface meeting the boundary, another interface (three or more regions) or a contact end. Inside a corner zone of radius 2s(x) around the feature, nodes are placed per sector (one per region, on the sector's bisector at distance s/2), and junction nodes on the boundary straddle the junction point. No angle limit is imposed. Instead the generator measures the dual interface near each feature (the union of Voronoi faces separating regions) and refuses (`degenerate_mesh`, feature index) when its Hausdorff distance from the true interface exceeds s, or a node falls in the wrong region. | Heuristic placement, checked postcondition (P4b). The bound s is the resolution M7 already accepts. | A check-based limit instead of an arbitrary angle. | Interface geometry inside corner zones is accurate only to s. |
| M9 | Contacts. 2D: a patch is a sub-interval of boundary segments; its ends are straddled along the boundary by node pairs (M8 junction rule). 3D: each convex face carries any number of named patch polygons (simple polygons in the face, pairwise disjoint, each inside one face), with the rest of the face under the face's name. The face is first meshed by the 2D generator, with its patch edges as straddled 1D interfaces; the face nodes then go into the 3D mesh. A node on an edge between two faces belongs to a patch of each, with each patch's own area. | Proven (needs M3 b): under boundary conformity a node's boundary area on a face is its 2D Voronoi cell among the face's nodes. So patch boundaries are exact to s, by the M7 argument in 2D. Gate: compare with the 2D path computed independently. | Planar MOSFET source, drain and gate on one top face. | A patch cannot span two faces; a contact over two faces lists both patches' nodes. |
| M10 | Slivers. **Not assumed harmless.** Facts: (i) exactly flat elements cannot exist (M4, P0); (ii) for Delaunay meshes the couplings and volumes are Voronoi measures, independent of slivers, and M1's exact determinants keep a sliver's circumcentre accurate; (iii) anything that interpolates on simplices (rendering, future AMR transfer, a future unstructured BTBT locator) **is** degraded by slivers. Handling: report the minimum dihedral angle and the sliver count (dihedral < 5°) per mesh. P5 compares the dual against an independent Voronoi clipping on test meshes. A sliver-rich stress test (perturbed lattice) must pass the same solver gates as a clean mesh. No refusal on angle. Sliver removal is deferred: it needs weighted Delaunay, so a power-diagram dual, a different control-volume definition needing its own ruling. | (i) and (ii) proven; their effect on the solve is measured, not assumed; (iii) limitation. | Separates what slivers can and cannot affect, with tests. | Simplex-based consumers must handle slivers when they arrive. |
| M11 | Refinement and termination. 2D: Ruppert/Shewchuk with concentric-shell segment splitting; minimum-angle bound default 20.7°, option up to 33.8°; skinny triangles allowed only at small input angles. 3D: radius–edge bound default 2.0 with encroachment splitting; option down to 1.414. A size field h(x) with h_min > 0 (else `invalid_input`) and grading ≤ 1.5 per h. Every run has a node cap (`MeshOptions::max_nodes`, default 2e6), giving `resource_exhausted` if hit. Modified steps (no Steiner point in a protected band: split the interface pair set locally instead; elements made only of pair nodes exempt from the quality bound) are heuristic. | 2D: termination at 20.7° for input angles ≥ 60°, cited and checked (Shewchuk 2002); below 60° concentric shells alone do not guarantee it (corrected 2026-10-15, U6), so such inputs are heuristic and capped. 33.8°: heuristic. 3D: radius–edge ≤ 2 for input angles ≥ 90°, cited and checked (Shewchuk 1998, U6); smaller input angles heuristic and capped. 1.414: heuristic (TetGen's practice). Band and pair modifications: heuristic, capped. | Proven where the literature proves it; otherwise capped and reported. | A heuristic run may stop at the cap, with a message naming it. |
| M12 | Postconditions, checked on every generated mesh in Debug and Release (O(n)): P0 every element has orient > 0. P1 locally Delaunay (exact, tie-broken). P2 every boundary segment or facet is Gabriel (exact). P3 M2 totals. P4a every edge with nonzero coupling joining two regions is a mirror pair whose midpoint lies on the interface within 1e-12·L_e (outside corner zones). P4b the M8 corner bound. P5 per-node closure: Σ A_e t̂_e (outward) + Σ boundary area · n = 0 within 1e-11 relative to Σ\|A_e\|. A failure is `degenerate_mesh` naming the check; it signals a generator bug, never a user error. | Checks, not claims. The 1e-11 is a few ulps times the face count (≤ ~200). | Each guarantee in M3, M7 and M8 is verified, not trusted. | ~10–20% on top of generation time (to measure). |
| M13 | Determinism and speed. Insertion order: BRIO rounds with a Hilbert-curve sort, using a fixed-seed splitmix64 (not `std::` distributions, whose output is implementation-defined). Point location by visibility walk from the last insertion (terminates on Delaunay triangulations); cavities by breadth-first search over neighbours. Same input and options give a bit-identical mesh on any machine (6.8). Budgets, measured and recorded, not promised: Release 2D 1e5 nodes ≤ 2 s; 3D 1e5 nodes ≤ 15 s. ctest meshes ≤ 1e4 nodes (2D) and ≤ 2e3 (3D); 3D solves in ctest ≤ 3e3 nodes; larger ones `[.benchmark]`. The Debug suite grows by ≤ 120 s in total, or cases move to `[.benchmark]`. | Expected O(n log n) for BRIO plus walking: cited (Amenta–Choi–Rote 2003). Budgets are targets. | 3D without fast location is O(n²); SparseLU fill-in makes large 3D solves slow (PARDISO for benchmarks). | Hilbert sort plus BRIO adds code. |
| M14 | Mesh metadata persists. The io mesh file gains optional keyed records (IO2: old readers skip them): simplices, node region, patch provenance, generator options and a generation digest. `mesh_digest` stays the graph only, so a mesh's identity is unchanged. A new `load_generated_mesh` returns `GeneratedMesh`, or `invalid_input` if the records are missing; `load_mesh` is unchanged. | — | Reloaded meshes stay renderable (Unit 25) and refinable (Unit 17) without regeneration. | Larger files (simplices are about 2× the edges). |
| M15 | (Ruled R4.) Doping profiles: constant, box, Gaussian and erfc per axis with lateral erfc fall-off, summed per species, sampled at nodes (`device/doping_profiles.hpp`). | — | Point sampling matches today's lambdas. | Abrupt doping steps are not straddled (U9). |
| M16 | Nonlocal BTBT stays tensor-only; it is already refused without `DeviceDescription::cells`. | — | Keeps the unit to meshing. | Unstructured nonlocal tunnelling waits for a cell-walk locator. |

**Implementation notes, commit 3 (predicates and simplex import; differences from the text above, for review):**
- **M4, exact stage:** exact integers over the inputs scaled by a common power of two, instead of Shewchuk's expansion
  arithmetic.
  - Why: expansions can underflow on the degree-6 diametral-sphere test and on circumcentre numerators at extreme
    coordinate scales. Integers have no range restriction and are simpler to verify.
  - The filter is a forward error bound computed alongside the same formula, instead of Shewchuk's stage-A constants, so
    the filter and the exact stage cannot evaluate different formulas.
  - Trade-off: the exact stage is slower; it runs only on near-ties.
- **M4, tie-break:** the largest node id is perturbed most. Any fixed order is consistent; this one is documented and
  pinned by a test.
- **M2, scales:** the scale is max(L^(D−1), Σ|pieces|) for couplings, max(l^D, Σ|pieces|) for volumes and
  max(l^(D−1), Σ|pieces|) for patch areas, instead of L^(D−1) and l^D alone. The Σ|pieces| term covers cancellation
  among signed pieces near poorly shaped elements, where rounding scales with the pieces rather than with L.
  - The P5 closure tolerance is relative to Σ|pieces| for the same reason.
  - A patch node whose area is within τ of zero is dropped from the patch, as zero couplings are dropped.
- **Import scope:** P0, P1, P2, P3 and P5 are checked, plus manifoldness and folds. Accepted, but stated as limits:
  - global overlap of a locally valid mesh is not detected;
  - P4 cannot be checked without interfaces (U2).
- **Determinism:** claimed only within one build. Topology decisions are exact and therefore platform-independent;
  bit identity of the measures across compilers or machines is not claimed, because it has not been verified.
- **Mutation tests:** 12 run against the 20 new test cases; 10 caught, 2 survived.
  - Caught: insphere sign; tie-break order; the query-point cofactor sign; the face-side sign of the 3D dual; the drop
    rule using `< 0`; P1 skipped; P2 skipped; the 3D boundary piece and the volume factor.
  - The filter with no error bound is caught too, in a second form: its first form failed to build on an
    unused-variable warning and is not counted.
  - **Survived (closed in round 2, below):** the product rule written as max + 1, and the underflow guard removed. The worst-case error they
    under-bound is not reached by any input the tests generate, so the filter's soundness rests on the error analysis
    in `predicates.hpp`, not on a test (U11).

**Approval blockers, round 2 (owner request 2026-10-14; commits `ce7e1cd`, `93a6d2b`, `dd32cbd`):**

*U11, the filter (closed).* The soundness argument is now written in full in `src/mesh/predicates.hpp`: the
arithmetic model (binary64, round to nearest, gradual underflow, no contraction), the induction for sums and
products, cancellation (the bound is relative to M, not to |v|), overflow (an infinity or NaN reaches v or m) and
underflow.
- **Defect found and fixed:** the filter flagged a product only when its magnitude underflowed. A product whose
  *value* underflows while its magnitude does not escaped the flag. This happens when a factor's value is small by
  cancellation, and it falls outside the rounding model the bound assumes. Both are now checked.
- **Removed as redundant:** a deferral on a subnormal bound. Every formula's root is a sum of products, and each
  product's m is 0 or at least DBL_MIN, so a nonzero bound is at least d·2^-1073. Its rounding error is then at
  most 1/(4d) relative, which the factor 2 covers.
- **Tests (`tests/mesh/filter_test.cpp`)** use a probe that exposes each evaluation's value, magnitude, depth, flag,
  bound and decision beside the exact value. They check:
  - the depth table;
  - |value − exact| ≤ bound decided in exact integers, for 600 near-degenerate samples per predicate at scales
    2^±40, with the threshold region (|v| within 4× of the bound) exercised for every predicate;
  - both underflow contracts, overflow, cancellation at offsets up to 2^52, subnormal bounds, and the arithmetic
    model itself.
- **Mutations:**

| Mutation | Result |
|---|---|
| U11-1 product depth max + 1 (the earlier survivor) | caught by the depth table |
| U11-2 underflow flag ignored (the earlier survivor) | caught |
| U11-3 value-underflow check removed (the new check) | caught |
| U11-4 magnitude-underflow check removed | caught, after a test for a product with an exactly zero factor value |
| U11-5 subnormal-bound deferral removed | redundant (argument above); the code was removed |
| U11-6 bound factor 2 → 0.25 | caught |
| U11-7 magnitude of a difference subtracts | caught |

*M4, the tie-break (verified).* The scheme documented in `predicates.hpp` is checked against itself evaluated
directly, independently of the implementation's cofactor expansion:
- lifted rows (p, |p|² + ε^(r+1), 1), with r the rank by decreasing id and ε = 2^-128, scaled to integers;
- their Leibniz determinant (24 or 120 terms), with incircle its sign and insphere its negation;
- on tied integer points on x² + y² = 625 and x² + y² + z² = 9, translated, with half the cases off the tie.
- All 600 2D and 300 3D cases agree (more than 100 and 50 of them exact ties).
- Mutations caught: tie-break order reversed, query cofactor sign, incircle query cofactor as a unit column.
- Guarantees:
  - every sign is exact for any finite doubles;
  - the perturbed predicates never return 0 and define one consistent regular triangulation;
  - on any platform satisfying the arithmetic model the same inputs give the same signs, so the same topology. This
    is by the argument, not by cross-platform runs (none were made).
  - The dual measures are floating-point; their bit identity across compilers or machines is not claimed.

*M2, the tolerance (verified).* `src/mesh/dual_tolerance.hpp` holds the whole rule. `make_dual_mesh` calls it for
every coupling, patch area and volume, always on the assembled total and its scale, never on a piece. Non-finite
totals are refused. Tests:
- exact boundaries with nextafter on both sides of ±τ·scale, at five scales, plus NaN and infinity;
- a near-cocircular square whose diagonal crosses the threshold at ε = 2e-12, under coordinate scalings from 1e-7 to
  1e5;
- Kuhn grids scaled from 1e-9 to 1e3 dropping the same edges;
- exact cancellation of ±4 pieces (dropped, not refused), and the exact P1 refusal that precedes M2 when the same
  quad is moved by one ulp;
- a −4 piece with a +6.68 piece, kept at its exact total;
- a kite whose pieces are 25× its edge, where only the Σ|pieces| term drops a resolved total of 1e-9.
- Mutations: 9 of 9 caught (both inclusive ends, non-finite, the volume bound, τ ×10 and ÷10, either term of the
  scale, the coupling scale's power).

*M7, the interface rule P4 (specified and implemented as a check; no generator).* `include/NiTCAD/mesh/interfaces.hpp`
gives the metric and its computation:
- deviation(F) = min over the declared pieces between the two regions of max over the face's vertices of the
  distance to the piece;
- the face's vertices are the circumcentres around the edge, and on the boundary the boundary-facet circumcentres
  and the midpoint;
- because F lies in the hull of its vertices and the distance to a convex piece is convex, this is an upper bound on
  the face's distance from the interface.

Acceptance:
- **P4a**, outside the corner zones (2 s around piece boundaries): 1e-12 L + 1e-13 X;
- **P4b**, inside the zones: s;
- every region pair that meets must be declared.

`check_interfaces` is the only way to establish P4. `make_dual_mesh` never does, and an imported multi-region mesh
fails it unless its nodes follow the rule (tested). Results on brute-force Delaunay layouts:

| Layout | Result |
|---|---|
| Flat, uniform or graded pairs; thin layer wider than both bands; 3D plane | exact (flat deviation 0) |
| Node inside the protected band; layer inside a band with staggered pairs; 3D broken mirror | refused |
| Node just outside the band (1.05 b) | accepted, as the band argument says |
| Right-angle corner, pairs 8 s away, no corner nodes | exact: by mirror symmetry the corner point is equidistant from both arms |
| Right-angle corner, bisector nodes, pairs graded to 1.5 s | 0.80 s |
| 60° corner, bisector nodes, graded pairs, oblique boundary junction | 0.95 s |
| T-junction | 0.31 s |
| Corner threshold | the 0.80 s layout is accepted at resolution 0.026 and refused by P4b at 0.024 |

M8 findings, which change what the generator must do:
- the bisector nodes make a right-angle corner worse, so the generator should not place them there;
- an oblique boundary junction needs the pairs *and* the boundary nodes graded toward it (1.5 s, 3 s). Without that
  the import failed P2, or the junction face deviated beyond s;
- below 60° nothing is measured, and 60° itself passes with a thin margin (U12).

Mutations: 12 run; 11 caught. M7-12 ("zero-coupling cross edges judged too") is **equivalent**: such a face is
degenerate, and its single point is the circumcentre of cells whose kept cross faces already contain it; if a cell
has no kept cross face, that point carries no interface measure.

**Owner rulings on U2 and U12 (2026-10-14):**
- **U2, P4 certification: required.**
  - A multi-region device (two or more regions) on a mesh that is not a tensor product must carry a valid P4
    certificate. A missing certificate is refused, and so is a stale one.
  - The certificate comes from `mesh::certify_interfaces`. It records digests of the dual mesh, the node regions and the
    declared interfaces; it is stale when the device's mesh or node regions differ from those it was issued for.
  - Tensor-grid meshes are exempt by provenance (owner's choice of the recommended option). `make_tensor_grid` tags its
    mesh, `Mesh::from_parts` verifies that a mesh tagged so has only axis-aligned edges, and io saves and loads the tag.
    A file written before the tag existed loads as a general graph.
  - Single-region devices need no certificate.
  - Enforced in `Device::create`, so every solve of such a device has one.
- **U12, corner rules: as proposed, with acute corners refused.**
  - A bend of exactly 90° between two interface segments gets no corner nodes, since mirror symmetry makes it exact.
  - Every other interface feature with all sectors at least 60° (bends, junctions of three or more regions, junctions
    with the outer boundary) gets graded pairs (1.5 s, 3 s) and one bisector node per sector at s / 2. It is accepted
    only if the generated mesh passes P4b; otherwise generation fails with `degenerate_mesh`.
  - A sector below 60° is refused (`invalid_input`) until measurements justify it.
  - The angle tests are exact (integer arithmetic on the input coordinates).

**2D generator, round 3 (2026-10-14; commits `740beba` to `8be2cc3`):**

*What was built:*
- **P4 certification (U2):** `Device::create` refuses a multi-region device on a 2D or 3D non-tensor mesh without a
  current certificate. Tensor grids are exempt by provenance (`MeshOrigin`, verified, saved by io). 1D meshes are exempt
  too: an interpretation of the ruling, for the same reason (every edge lies along x).
- **The generator** (`mesh::generate_mesh`):
  - exact validation and the U12 corner rules;
  - mirrored pairs graded toward features, a straddle at boundary junctions, nodes exactly at patch ends;
  - Bowyer–Watson with the exact predicates, then Ruppert refinement with protected bands;
  - import and certification of the result.
- **Doping profiles (M15) and generated meshes on disk (M14).**

*Measured:*

| Gate | Result |
|---|---|
| Generator validity, quality (min angle ≥ 20.7° outside bands), size, determinism, refusals | passed |
| P4 by feature: flat, 90° bend, T-junction, 60.3°, 120°, thin layer, band-forcing refinement | certified: corner deviation 0, 0, 0.31 s, 0.95 s, 0.73 s |
| Speed (Release) | 103,627 nodes in 1.70 s (target 2 s) |
| Linear-field consistency (oxide between electrodes, skew mesh) | linear to 1e-10 V: passed |
| Smooth diode convergence | H1 1.79, 1.24; L2 2.36, 2.11; current 1.50, 2.09: passed (≥ 0.9, 1.5, 0.9) |
| Abrupt diode against 1D | error 1.9e-4 (11,549 nodes) against the tensor grid's 2.6e-4 (11,536): passed |
| MOSFET, k = 1.5 against the tensor MOSFET | **not met**: on current 4 % low at 8,031 nodes, tensor 0.05 % at 9,490 (U14) |

*Mutations:*
- certification: 7 of 7 caught;
- generator: 15 of 15 caught (one rule, boundary grading at junctions, was removed as redundant after its mutant
  survived: refinement splits do that work);
- doping and io: 7 of 7 caught.

*Findings that need a ruling:*
- **U14, the MOSFET gate (unmet).** The generator makes isotropic triangles. The tensor MOSFET resolves the inversion
  layer with 2e-8 cm cells normal to the surface; isotropic triangles that fine over the channel exceed the node budget.
  - Options: anisotropic refinement near selected boundaries (a boundary-layer mode for the generator); a hybrid mesh
    with a tensor strip under the gate; or a ruling that changes the gate.
  - The gate stays in the suite as `[!shouldfail]` with its measurement.
- **U15, the U12 corner rule.** Bisector nodes make every measured corner worse:

  | Corner | With bisector nodes | Without |
  |---|---|---|
  | T-junction | 0.31 s | 0 |
  | 60.3° | 0.95 s | 0.36 s |
  | 120° | 0.73 s | 0.32 s |

  Proposed: drop the bisector nodes everywhere. Graded pairs alone give a smaller deviation and a wider margin at 60°.
  They are implemented as ruled until then.
- **U16, U2 in 1D.** The 1D exemption above was applied without an explicit ruling.

**Owner rulings on U14, U15 and U16 (2026-10-14), and their implementation (commits `d741912`, `898fbf1`, `13f9148`):**

- **U15, corners (ruled: remove the bisector nodes where they worsen accuracy). Done.** They worsened every measured
  configuration, so none are placed. P4b still judges every feature, the 90° rule stands, and sectors below 60° are
  still refused. Regression tests pin each measured configuration, with deviation bounds from the measurements and no
  node within s of the feature:

  | Configuration | Corner deviation |
  |---|---|
  | Flat interface | 0 |
  | 90° bend | 0 |
  | T-junction | 0 |
  | 60.3° | ≤ 0.37 s (measured 0.359) |
  | 120° | ≤ 0.33 s (measured 0.321) |

  Mutations: 3 of 3 caught, including the old bisector generator run as a mutant.
- **U16, 1D (ruled: exempt only verified 1D tensor grids). Done.**
  - The blanket 1D exemption is removed.
  - The tensor tag is verified by `Mesh::from_parts`; in 1D it also requires a chain: strictly increasing x, and
    exactly the n − 1 edges joining x-consecutive nodes. Every region transition then lies on the dual face between
    two neighbours.
  - Any other 1D mesh with two regions needs `mesh::certify_interfaces_1d`: every cross-region edge's midpoint (its
    dual face) on a declared interface point of that pair, within the P4a tolerance.
  - The 1D insulator test, a mesh built from parts, now carries such a certificate.
  - Mutations: 5 of 6 caught. The sixth, the increasing-x check removed, is equivalent: with n − 1 rank-consecutive
    edges, a duplicate x is an edge of zero length, which the length check refuses.
- **U14, the MOSFET gate (ruled: keep k = 1.5; investigate anisotropic layers and a hybrid tensor strip; pass before
  3D). Passed.**
  - Both options were realised as one mechanism, `MeshOptions2D::layers`: a structured boundary layer is a tensor
    lattice along a boundary segment, and placed under the gate it is the hybrid strip. Its cells are Delaunay ties,
    so its dual is exactly the tensor box method; refinement stays out of it.
  - Mesh specifications tried, the threshold fixed throughout:

    | Attempt | Specification | Nodes | Result |
    |---|---|---|---|
    | 0 | isotropic, size boxes | 8,031 | on current 4 % low |
    | 1 | layer to 1e-6 cm | 4,438 | subthreshold 1–4 % high |
    | 2 | layer to 1.2e-5 cm | 7,337 | high-gate current and swing fail |
    | 3 | layer to 1.6e-5 cm | 8,304 | high-gate current, g_m and swing fail |
    | 4 | as 3, contact edges exactly at Lsd and Lsd + Lg in all three meshes | 8,190 | only the swing fails, 7.5e-4 against 1.03e-4 allowed |
    | 5 | layer to 1.8e-5 cm | 8,864 | swing fails, 1.4e-4 |
    | 6 | layer to 1.8e-5 cm, unstructured refinement (1e-6 cm) below at the body contact | 9,342 | **passes** |

  - Before attempt 4, the tensor baseline and reference missed the contact edges by 2.9e-8 and 1.4e-8 cm, an offset
    that halves with refinement. Snapping their nearest nodes onto the edges, with node counts unchanged, makes all
    three meshes model the defined device. It does not weaken the rule: the tensor error only shrinks.
  - **Plainly:** the passing mesh is a tensor lattice over 90 % of the depth with unstructured refinement only at the
    bottom. Isotropic triangles alone could not meet the gate at this node budget.
  - The gate is in ctest (`solve_unstructured`, Release) as an ordinary test; the `[!shouldfail]` marker is removed.
  - Mutations: 8 of 9 caught. The ninth, judging layer-only triangles, is equivalent: they are cell halves whose
    circumcentres lie in the layer, already exempt.
  - A defect found and fixed: boundary nodes were filled between lattice columns more than h apart, breaking the
    lattice. Hidden in the MOSFET (columns closer than h), found by the new lattice test.

**3D generator, commit 6 (owner approval 2026-10-15; commits `c12915a`, `86621b7`, `9da439a`, `3c56d24` and the 6e
gates):** Owner choices for this commit: full M6 (convex polyhedra), planar interfaces with no junction lines, and a
z-uniform convergence gate against a 2D reference.

- **6a, `Triangulation3D`** (`src/mesh/triangulation3d.hpp`):
  - Bowyer–Watson insertion with the exact orient3d and the perturbed insphere, and a visibility walk.
  - BRIO insertion order: splitmix64 rounds, with a 3D Hilbert order inside each round (Skilling).
  - It equals the brute-force perturbed Delaunay on a lattice where every tie occurs.
  - 1e5 points take 4.2 s (Release).
  - Defect found and fixed: a cavity slot reused before its outside neighbours were relinked.
- **6b, single region** (`generate_mesh(Geometry3D, MeshOptions3D)`):
  - Exact validation of convexity, closure and planar, convex faces; non-convex domains are refused (U7).
  - Refinement follows Shewchuk:
    - encroached domain-edge subsegments are split, with concentric shells at vertices;
    - encroached subfacets are split at their circumcentres;
    - tetrahedra whose radius–edge ratio exceeds the bound, or whose circumradius exceeds √(3/8)·h, get their
      circumcentre unless it encroaches.
  - **Change from the plan:** faces are not meshed by the 2D generator. Each face keeps a 2D Delaunay triangulation of
    its nodes, in face coordinates and exact on axis-aligned faces, and that triangulation defines its subfacets. The
    public 2D generator is unchanged, and the face and 3D nodes are one set.
  - Extraction:
    - drops the super-vertex tetrahedra and any tetrahedron whose four vertices lie on one face;
    - checks that every exposed triangle lies on one face and that the volume matches to 1e-12;
    - then runs `make_dual_mesh`.
  - **Oblique faces pass**, including an irregular tetrahedron at device scale whose split points round off its lines and
    planes. The pre-stated fallback (refuse oblique faces, U17) is not needed, and no flat tetrahedron occurred.
- **6c, face patches (M9)** (`PolyhedronFace::patches`):
  - A conforming tiling of a face. Its interior edges become constrained subsegments, so every boundary triangle lies in
    one patch and patch areas are exact.
  - **Change from M9:** M9 said patch edges were straddled and accurate to s. Conforming edges make them exact, as the
    2D generator already does at patch ends.
  - Gate: every node's patch area matches `make_dual_mesh` 2D on that patch's triangles to 1e-12. Cases:
    - contact strips;
    - a contact inside a face;
    - a slanted patch edge.
- **6d, planar interfaces (M7)** (`Geometry3D::interfaces` and `regions`):
  - **Validation:**
    - each plane stays at least 4 s from every vertex;
    - it meets every face it cuts at 60–120°, checked before its traces;
    - each trace meets the domain edges it crosses at 60–120° and is at least 6 s long;
    - no cut face carries patches (deferred);
    - planes never meet in the domain and stay at least 4 s apart;
    - there is one seed per slab.
  - **Nodes:**
    - on each crossed edge, nodes at ±s/2 from the crossing (the segment between them is never split);
    - along each trace in a face, pairs at ±s/2 within the face, graded from each end (1.5 s, 3 s);
    - on the polygon inset by 1.5 s, a pair set P from the 2D generator, with nodes at p ± (s/2)·n.
  - **Band:** b = √(h_P²/3 + s²/4). Inside it:
    - no circumcentre is inserted;
    - an encroached subfacet gets a new trace pair;
    - tetrahedra made only of interface nodes are exempt.
  - The result carries the P4 certificate. `Device::create` refuses the mesh when the certificate is missing or stale.
  - P4b measured at the traces (s = 0.02, h = 0.25, unit cube), with regression bounds taken from the measurements:

    | Configuration | Flat deviation | Corner deviation | Bound |
    |---|---|---|---|
    | Horizontal plane | 0 | 0 | 1e-12 s |
    | Thin middle layer (5 s) | 0 | 0 | 1e-12 s |
    | Oblique plane | 3.75e-3 of tolerance | 0.180 s | 0.19 s |

- **6e, gates** (predefined; none set after a run):
  - **Linear field:** conjugate gradients on the generated dual give the linear potential to 1e-12.
  - **Kuhn equivalence:** a solve on Kuhn-split boxes matches the tensor diode solve to 1e-12 in potential and current.
  - **io:** a 3D generated mesh round-trips with its patches, regions, certificate and the 3D report (optional records).
  - **Convergence** (z-uniform doping; reference a 2D tensor solve 4× finer than the finest mesh; `[.unstructured3d]`,
    beyond the ctest node budget of M13, run explicitly):

    | Mesh | Nodes | H1 orders | L2 orders | Current orders |
    |---|---|---|---|---|
    | Generated | 391, 2030, 13057 | 1.78, 1.32 | 1.63, 2.05 | 1.63, 2.18 |
    | Perturbed lattice (sliver stress, M10) | 132, 693, 4305 | 2.12, 1.81 | 2.62, 2.31 | 2.22, 2.01 |

    - The generated meshes have 6, 13 and 133 slivers, with minimum dihedral angles down to 0.65°. The perturbed
      lattices have 0, 33 and 290, down to 0.24°. Slivers did not degrade any gate, which M10 (ii) predicts and the run
      now measures.
    - One of my own test assertions failed: every stress mesh must have slivers. The coarsest lattice has one cell
      through z, so no interior point moves and it has none. I narrowed that check to meshes with interior points. The
      gates were unchanged.
- **Mutations** (each run against the tests; "equivalent" means it survived every measured case, with the reason):

  | Commit | Caught | Equivalent in the measured cases |
  |---|---|---|
  | 6a | 3 of 3 | — |
  | 6b | 2 of 4, plus a fault injection caught by the postcondition and, with the postcondition removed, by the tests | skip subsegment encroachment (subfacet splitting splits the subsegments it encroaches); remove the flat-tetrahedron drop (none occurs) |
  | 6c | 4 of 5 | skip patch-segment encroachment (same reason as 6b) |
  | 6d | 6 | band rule alone and pair-only exemption alone (they mask each other; caught together); trace pairs at ±3s/4 (still mirror-symmetric, so exact); trace-edge rule loosened (on these shapes the face rule refuses the same planes) |
  | 6e (io) | 2 of 2 | — |

- **Performance** (Release, budget targets, M13):
  - As first committed: 128,671 nodes took 23 s, the 15 s per 1e5 target missed. `make_dual_mesh` alone took 16.5 s,
    14.3 s of it in the exact coupling pieces (M1).
  - After U18 (below): 15.9 s, about 12.4 s per 1e5, so the target is met; `make_dual_mesh` takes 7.8 s.
  - The Debug mesh tests grew to 61 s; the 3D tests add about 33 s.
- **Literature (U6), checked 2026-10-15 (completed the same day):**

  | Source | What was checked | Status | Consequence here |
  |---|---|---|---|
  | Gärtner–Si–Fuhrmann 2010 (Comput. Math. Math. Phys. 50(1)) | Abstract | **Confirmed.** 3D boundary-conforming Delaunay meshes when the smallest input angle is at least arccos(1/3) ≈ 70.53°; arbitrary 2D polygonal regions | M3 (b), M5 and M6 cite it as stated |
  | Shewchuk 1998 (SoCG, pp. 86–95) | The paper itself. Its math glyphs are CCITT bitmaps; I decoded them, and the same glyphs read 90° in Ruppert's 2D condition in that paper, which the literature gives as 90° | **Confirmed.** Radius–edge ratio ≤ 2 when no input angle (segments and facets, dihedral included) is below 90° | 3D termination is cited only for domains whose input angles are all at least 90° (boxes). The tetrahedron (70.53°), prism (60°) and octahedron (60° face angles) cases, patch edges and interfaces are heuristic and capped |
  | Shewchuk 2002 (the 2001 report of the CGTA article) | The paper itself, glyphs decoded | **Confirmed:** Ruppert's algorithm terminates when no input angle is below 60° and the split threshold is B ≥ √2, so no angle ends below arcsin(1/(2√2)) ≈ 20.7°. **Corrected:** with smaller input angles, plain Ruppert refinement, concentric shells included, can fail to terminate; only the "Terminator" (shells plus a look-ahead) guarantees termination | M11's "with small angles under concentric shells" was wrong. The 2D generator has shells but no look-ahead, so inputs with angles below 60° (domain corners; interface sectors below 60° are already refused) are heuristic and capped by `max_nodes` |
  | Eymard–Gallouët–Herbin 2000 (Handbook of Numerical Analysis VII) | Through the same authors' 2024 paper (arXiv 2405.16985, introduction) | **Confirmed second-hand:** order-1 error for the two-point flux scheme on admissible meshes when the solution is in H². The handbook's discrete-H1 statement itself was not reached | The H1 ≥ 0.9 gate rests on an order-1 result |
  | Amenta–Choi–Rote 2003 (SoCG, pp. 211–219) | Abstract | **Confirmed:** BRIO is no worse asymptotically than a random insertion order | Visibility-walk location is outside that analysis; the measured 3D Delaunay speed (1e5 points in 4.2 s) is the evidence here |
  | Devillers–Teillaud 2011 (CGTA 44:160–168) | Bibliography only; the HAL copy was blocked | **Not reached** | The tie-break does not depend on the citation: tests check it against the documented lifted perturbation, brute-force Delaunay and antisymmetry |

**Gates (predefined, so a threshold is never set after the run):**
- **Tensor equivalence:** split rectangles and Kuhn-split boxes give exactly `make_tensor_grid`'s volumes, couplings
  and areas to 1e-13 relative, with the diagonals dropped. Their diode solve matches the tensor run to 1e-12.
- **Dual correctness:** P5 closure, plus an *independent* Voronoi cell per node (half-space clipping by bisectors, which
  does not use the simplices) equal to the simplex-built dual within 1e-12 on small 2D and 3D meshes. Closure alone
  cannot detect a wrong dual that is not orthogonal; the independent cells can.
- **Linear-field consistency at the assembler:** Poisson with no charge, Dirichlet at two opposite faces and Neumann
  elsewhere. The discrete potential equals the exact linear one to 1e-12 on generated meshes, and the field
  reconstruction recovers the uniform field exactly. This replaces revision 1's linear-flux gate, which duplicated
  closure.
- **Convergence** (smooth Gaussian doping, h, h/2, h/4; reference a 4× finer tensor solve interpolated
  bilinearly):
  - discrete H1 potential error order ≥ 0.9 (proven O(h) for admissible meshes: Eymard–Gallouët–Herbin, cited, to
    verify);
  - L2 potential order recorded, gate ≥ 1.5 (observed second order is typical but not proven on general meshes);
  - terminal current order ≥ 0.9.
  - Abrupt junction: a laterally uniform 2D unstructured diode against the 1D legacy-validated diode; current error ≤
    the error of a 2D tensor grid with the same node count against the same reference.
- **MOSFET acceptance:**
  - Reference: the existing tensor MOSFET (`tests/solve/mosfet_test.cpp`) at 2× refinement per direction.
  - Compared: I_D at every sweep point above the current resolution, V_T (analysis extraction), maximum g_m, and the
    subthreshold slope over the resolved range.
  - Pass: |Q_unstructured − Q_ref| ≤ k · |Q_tensor − Q_ref| + f, with k = 1.5 (U1). Q_tensor is the baseline tensor
    MOSFET; the unstructured mesh has at most as many nodes. f = the reference's Richardson error estimate.
  - Every point converges, with bias step cutting on if needed (recorded).
  - Off-state current is excluded (6.13).
- **Slivers:** the perturbed-lattice 3D diode passes the same convergence gates as the clean one; sliver statistics are
  recorded.
- **Refusals, predicates, determinism, io:** each M2/M6/M8/M12 refusal returns its code and index; exact signs on 1e-17
  perturbations and cocircular/cospherical sets; two runs bit-identical; `save`/`load_generated_mesh` round trip.
- **Mutation tests:** at least: drop a circumcentric piece's sign; skip P2; skip the band rule; reverse the tie-break; use
  `<` for `<=` in τ; break insertion order. Each must fail a gate.

**Unresolved (need a ruling or carried as limits):**
- U1 MOSFET factor k = 1.5: a judgement, not derived. The owner may set it.
- U2 **Ruled 2026-10-14** (above): multi-region devices on non-tensor meshes need a valid P4 certificate.
- U3 Termination is unproven for the band and pair modifications (M7, M8), for 3D convex domains with an
  input angle below 90° (Shewchuk 1998's hypothesis, checked in U6), for 2D inputs with an angle below 60° (Shewchuk 2002), and for options beyond the cited bounds. Capped (M11) and reported;
  every 3D test case terminated well below its cap.
- U4 Corner-zone accuracy (M8) has no a-priori bound; it is measured and refused above s (P4b). Measured in round 2: 0 at a
  right angle without corner nodes, 0.80 s with bisector nodes, 0.95 s at 60°, 0.31 s at a T-junction.
- U5 Slivers affect future simplex-based consumers (M10 iii); removal needs a power-diagram dual (its own decision).
- U6 **Checked 2026-10-15** (table in "3D generator, commit 6"). Confirmed: Gärtner–Si–Fuhrmann (70.53°), Shewchuk 1998
  (ratio ≤ 2 for input angles ≥ 90°), Shewchuk 2002 (20.7° for input angles ≥ 60°), Amenta–Choi–Rote (BRIO); second-hand:
  Eymard–Gallouët–Herbin (order 1, u in H²). Corrected: concentric shells alone do not guarantee termination below 60°
  (M11). Not reached: Devillers–Teillaud (the tie-break is tested independently).
- U7 Non-convex 3D domains (fins, mesas) are deferred. So are 3D interface junction lines and bends (owner choice
  2026-10-15), interfaces cutting faces that carry patches, and 3D structured boundary layers.
- U8 The L2 order gate (1.5) is a judgement where theory gives no rate.
- U9 Region-only straddling: doping steps (abrupt junctions not on a region boundary) are not straddled, so an abrupt
  junction is resolved only to h. "Doping interfaces" with pairs could be added later.
- U11 Filter soundness: **closed** in round 2 (the argument written out, one defect fixed, direct tests against the exact
  value, every mutation caught or removed as redundant).
- U12 **Ruled 2026-10-14** (above): no corner nodes at 90° bends; graded bisector nodes elsewhere, accepted only if P4b
  passes; sectors below 60° refused.
- U13 P4 does not check that each node lies inside its region's geometry, which is not an input. A misplaced node
  surrounded by its own region's nodes is not detected; one next to another region is.
- U14 **Resolved 2026-10-14**: the MOSFET gate passes with a structured boundary layer (the hybrid strip), k = 1.5
  unchanged (above).
- U15 **Resolved 2026-10-14**: no bisector nodes; regression tests per measured configuration (above).
- U16 **Resolved 2026-10-14**: only verified 1D tensor chains are exempt; other 1D meshes need certify_interfaces_1d.
- U10 Budgets (M13) are targets. If the 3D generator misses its gates or budgets within the unit, the PR carries 2D plus
  3D import and says so (unchanged from revision 1). Result: every 3D gate passed; the 3D speed target was
  missed at first (23 s at 1.29e5 nodes) and met after U18 (15.9 s).
- U17 (pre-stated fallback for oblique faces) **not needed**: oblique faces pass.
- U18 **Resolved 2026-10-15** (owner request), with no change to M1.
  - Profile (128,671 nodes, 774,877 tetrahedra): the exact coupling pieces took 14.3 s of `make_dual_mesh`'s 16.5 s,
    the exact circumcentres 3.1 s, the predicates 0.05 s. The cost was the exact integers' storage, not their
    arithmetic: each operation allocated its limbs on the heap.
  - Fix: `ExactInteger` keeps up to 24 limbs inline and moves to the heap only beyond that; its loops use raw
    pointers; `to_double` gathers its top 64 bits from at most three limbs, with the same truncation as before.
  - **Bit-identical:** the same input gives the same bits. Over all 774,877 tetrahedra, the digests of every
    coupling piece and circumcentre, and the dual graph digest, are equal before and after the change.
  - Result: couplings 14.1 s to 5.6 s; `make_dual_mesh` 16.5 s to 7.8 s; the 3D benchmark 23 s to 15.9 s at
    1.29e5 nodes (target met); the 2D benchmark 1.96 s to 1.15 s at 1.04e5 nodes.
  - Tests: inline and heap values to about 94 limbs; `to_double` against a known top 64 bits at every shift
    offset; identities through heap-sized products.
  - Mutations: 4 of 5 caught. The fifth, heap growth not copying the limbs, is equivalent: growth only ever
    copies limbs that `assign` has just zeroed, into a zeroed buffer.
  - The earlier idea, a filtered circumcentre certifying the correctly rounded value, is not needed; it would have
    changed M1, since `to_double` is faithful, not correctly rounded.

