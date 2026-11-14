<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 25c — vector glyphs, the clipping box, the exploded view, transient filmstrips (2026-10-22; on `render/vectors-clip`)

**Owner request:** VQ9 (2026-10-17): "Unit 25c (vector glyphs, clipping box, exploded view, transient filmstrips)
after 25b, each view droppable"; the plan's views 20 (current arrows on a 2D section), 21 (E arrows on the surface),
23 (a clipping box, a quarter cut away), 24 (regions and electrodes lifted apart) and 26 (a filmstrip of n as the
gate switches on), about 1,600 lines. On 2026-10-22 the owner asked for 25c to start without a further approval
("then start 25c without my approval"): the rulings below (CQ1–CQ4) are taken as recommended under that instruction
and recorded for the owner's review at acceptance.

**What exists (audited on `main` at `8256cf8`).** Nodal E and J and the no-flux streamline field (Unit 25b), point
location, slices and convex-polygon clipping (`Geometry::slice`, the internal `clip`), the 3D solid layer by
control-volume pieces with translucent insulators, contact slabs, the outline shell. Transient results hold
snapshots with fields, edge currents and rates, but `nitcad plot` converts a snapshot to a bias point with its fields
and bands only, so \|J\|, streamlines and R are unavailable on transient results (finding; fixed here).

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| DC1 | **Vector glyphs** of E, J, Jn or Jp (`vectors.hpp`'s streamline field, so tangential on no-flux boundaries): 2D at the points of a regular grid over the device, about N across; 3D at points of the boundary surface, spread evenly by area (a deterministic stratified sample of its faces), or on a slice. Each arrow a flat shaft and head (2D, over the map) or a lit cone on a cylinder (3D), of one length (a fraction of the spacing), colored by \|v\| through the colormap on the field's scale; points with no field drawn nothing. | Direction is what arrows add over a map; one length keeps 20 decades of magnitude readable, and the color carries the magnitude (CQ1). | No length-coded magnitude. |
| DC2 | **Clipping box** (3D): an axis-aligned box removed from the device. The remaining boundary faces (each clipped to its parts outside the box, as convex pieces) and the cut faces (the box's faces inside the device, from slices of the geometry clipped to the face) are drawn with the field or the regions, insulators translucent as in the field view. | The standard way to see inside a 3D device while keeping its outline; exact from the existing convex clipping. | The box is axis-aligned. |
| DC3 | **Exploded view** (3D): each region's closed surface (its boundary faces and its side of every interface) and the contact slabs on it, moved apart along an axis (default y, the stack) by a gap per step in the regions' order along it; drawn solid by region or with a field. | Shows the layer stack and buried regions without slicing. | Regions that interlock along the axis still overlap in projection. |
| DC4 | **Transient filmstrip**: chosen snapshots (all, or every k-th, at most 12) drawn as one field view each and tiled into one image, each titled with its time, on one scale fitted over all of them (CQ2); `--frames <dir>` writes each frame as its own PNG for a video tool. | Comparing frames needs one scale; frames serve presentations. | No animation playback (26). |
| DC5 | **`nitcad plot`**: `--vectors E\|J\|Jn\|Jp [--arrows n]` over any 2D or 3D field or slice view; `--clip <x0>,<y0>,<z0>:<x1>,<y1>,<z1>` (µm) on 3D field and solid views; `--view exploded [--gap <um>] [--axis x\|y\|z]`; `--view filmstrip [--every k] [--columns c] [--frames <dir>]` on a transient result. Transient snapshots keep their edge currents and rates. | One flag per view on the existing command. | — |
| DC6 | **Out of scope:** interactive clip handles, animation (26); oblique clipping planes; streamtubes; glyphs scaled by magnitude; video encoding. | Keeps 25c to VQ9's list. | — |

**Rulings (CQ1–CQ4), taken as recommended under the owner's instruction of 2026-10-22:**

| # | Question | Ruling |
|---|---|---|
| CQ1 | Arrow length by magnitude, or one length with the magnitude as color? | **One length, colored by \|v\|** on the field's scale. |
| CQ2 | Filmstrip frames on one shared scale, or each fitted? | **One scale** fitted over all chosen frames. |
| CQ3 | Exploded gap default? | **12% of the device's extent along the axis**, per region step. |
| CQ4 | Clipping applied to which views? | **Solid and field views in 3D** (and the vectors over them); slices and isosurfaces unchanged. |

**Acceptance gates:**

| # | Gate |
|---|---|
| CG1 | **Glyphs:** in a uniform field every arrow is parallel to it and of the set length; 2D arrows lie on the grid's points inside the device, as many as asked (to ±10%); 3D arrows lie on the boundary and are tangential on no-flux faces; colors are the colormap at \|v\|'s position; a NaN region gets none. |
| CG2 | **Clipping:** for a box inside a tensor device, the kept boundary pieces and the cut faces together form a closed, outward surface whose enclosed volume is the device's minus the box's (to rounding); a box covering a corner (the quarter cut) likewise; a box missing the device leaves the view unchanged; field values on the cut faces equal the interpolated field. |
| CG3 | **Exploded:** each region's surface is closed and outward; its volume is the region's control volume total; each region and its contacts move by the stated offset; with a zero gap the surfaces cover the solid view's. |
| CG4 | **Filmstrip:** one tile per chosen snapshot in time order, one shared scale, titles with the times; `--frames` writes one PNG per frame; transient snapshots draw \|J\| and streamlines. |
| CG5 | **Images:** WARP references (25a's tolerance): the 2D MOSFET's \|J\| with current arrows, the 3D MOSFET's E arrows on its surface, a quarter-cut clipping box of n, the exploded MOSFET, a filmstrip; a clean D3D12 debug layer in Debug. |
| CG6 | **`nitcad plot`:** each new flag from decks and result files; refusals with a message (exit 1, or usage 64): vectors in 1D, a clip box in 2D or not a box, a filmstrip of a non-transient result, unknown axes. |
| CG7 | **No numerical change** (only `render`, `nitcad plot`, tests and documents change); Release and Debug suites at every commit, no warnings under `/W4 /WX`; mutants in the arrow placement and direction, the clipping (inside/outside, the cut faces), the explosion offsets and the filmstrip scale and order, all caught. |
| CG8 | **Benchmarks:** arrows, clipping and the exploded view on the 10⁶-tetrahedron device, recorded. |
| CG9 | **No scope creep:** nothing outside DC1–DC5; ARCHITECTURE.md 6.14 and 11.3 as built after acceptance. |

**Commits (planned):** 1) vector glyphs and `--vectors` (CG1); 2) the clipping box and `--clip` (CG2); 3) the exploded
view (CG3); 4) filmstrips, transient snapshots' currents and rates in the plot, references and benchmarks (CG4, CG5,
CG8). Each view is independent and can be dropped (VQ9).

**Commit 1 (vector glyphs; CG1, CG6 in part), as built.** `streamlines.hpp`: `Arrow` (site, unit direction, \|v\|,
the surface's normal, the interpolant), `grid_arrows` (2D: cell centres of a square grid centred on the mesh's box,
spacing sqrt(area / count)), `surface_arrows` (3D: a stratified sample of surface polygons by area, each with its
polygon's normal); both read the field at the site through the point locator and skip undefined or zero fields.
`scene.hpp`: `add_arrows` (2D a flat shaft and head just in front of the map, unlit; 3D a six-sided shaft and cone
lifted a quarter of the length off the surface, lit; one length; one color or \|v\| through a colorbar),
`exposed_boundary` (the boundary faces not wholly on a contact). `nitcad plot --vectors E|J|Jn|Jp [--arrows n]`:
arrows of the streamline field (Unit 25b: tangential on no-flux boundaries), 0.8 of the spacing long; 2D over the
map (by default the vectors' own magnitude) in near-black; 3D, over their own magnitude, on the solid device colored
by \|v\| with its colorbar (CQ1), else near-black over the chosen field or slice view; 3D sites on the exposed
boundary or the slices shown. **CQ1 as built:** colored by \|v\| where the view does not itself show a field (the
solid); over a field map an arrow colored like the map under it would vanish, so there it is near-black.
- **CG1:** a uniform field's grid arrows (50, 200 and 800 asked) number within 10% of the request, point along the
  field (10⁻¹²), lie on a grid of the stated spacing centred on the device; none where the field is NaN or zero;
  surface arrows on a 3D bar's exposed faces (its contact faces excluded: 6 sides less 2) have unit normals and no
  normal component away from the contacts; drawn arrows span the set length along their direction, sit in front of
  the 2D map or above their 3D surface, and carry the color or \|v\|'s colormap position.
- **CG6 (in part):** J arrows in 2D and 3D, E arrows over n (2D) and over a slice (3D), Jn arrows with streamlines;
  refused: an unknown field and `--arrows 0` (usage), a 1D device.
- Looked at on a 2D and a 3D MOSFET deck (the 3D one at V_G = 5 V, V_D = 0.5 V): current into the drain pocket and
  out of the source, the body leakage, fields from the gate; arrows on contact faces stood normal to them under the
  slabs, so contact faces are excluded (`exposed_boundary`).
- **Mutation checks (CG7, Release): 11 of 11 caught** — zero fields kept, the direction unnormalized, the spacing
  doubled, the grid not centred (caught once the centring check was added), a non-unit normal, the 2D arrow behind
  the map, a short 2D and 3D arrow, the 3D arrow sunk into its surface, contact faces kept, the color not \|v\|'s.

**Commit 2 (the clipping box; CG2, CG6 in part), as built.** `geometry.hpp`: `Box`; `outside_box` (a convex
polygon cut by the box's six planes in turn, the part beyond each kept and the rest passed on; a part lying in a box
plane, with nothing strictly beyond it, belongs to the box); `Geometry::clipped_surface` (the boundary faces' parts
outside the box, and the cut faces: sections of the elements by each box face strictly inside the mesh's bounding
box, clipped to the face and turned toward the box, each with its interpolants). `ViewOptions::clip` (CQ4: the 3D
solid and field views draw the clipped surface, interfaces cut alike, a contact slab's face left out when its
centroid is in the box). `nitcad plot --clip <x0>,<y0>,<z0>:<x1>,<y1>,<z1>` (µm), arrows then sit on the clipped
view's exposed boundary.
- **CG2:** a square's parts outside a box are convex and cover the square less the overlap; on a graded tensor mesh,
  an interior box, a corner box (the quarter cut), a slab through the device (two pieces), a box with a face on the
  boundary and a box missing the device give closed surfaces (vector area 0 to 10⁻¹⁰) enclosing the device's volume
  less the overlap (10⁻¹⁰), outward; the same for an interior and a corner box on a generated mesh; a linear field
  exact on the cut faces; the field and solid views (with oxide above and below, so interfaces too) draw nothing
  inside the box and draw the cut face; a box over the drain end removes its pad.
- **CG6 (in part):** clipped field, solid and arrow views from a deck, and the clipped image differs from the whole;
  refused: a two-coordinate or inverted box (usage), a slice view, a 2D device.
- Looked at on the 3D MOSFET deck: the quarter cut shows n inside the source and drain pockets and the channel
  region, the gate slab notched where the box removes it.
- **Mutation checks (CG7, Release): 10 of 10 caught** — the box's faces swapped, in-plane parts kept, the remainder
  not reduced, cut faces at the bounding box, cut faces turned outward, cut faces not clipped to the face, slabs
  kept, the view unclipped, interfaces unclipped, `nitcad plot` ignoring the box (caught once the clipped and whole
  images were compared). The face-on-boundary case was added after a first version kept such a face.

**Commit 3 (the exploded view; CG3, CG6 in part), as built.** `geometry.hpp`: `Geometry::region_surfaces` (per
region, the boundary pieces its nodes own and its side of every interface piece, each turned outward from the
region). `scene.hpp`: `exploded_view` (each region's surface and its contacts' slabs moved along the axis by the
gap times the region's rank in the order of its mean node coordinate, the first staying; the semiconductor with a
field or by class, insulators translucent), `default_gap` (CQ3: 12% of the extent). `nitcad plot --view exploded
[--field f] [--axis x|y|z] [--gap um]` (regions without `--field`).
- **CG3:** on a 3D silicon device under an oxide layer with an electrode, each region's surface is closed (vector
  area 0 to 10⁻¹⁰) and outward, enclosing exactly its nodes' control volumes (10⁻¹⁰); at zero gap the surfaces' area
  is the boundary's plus both sides of the interface; the oxide and its electrode move up by one gap and the silicon
  stays; a field colors the silicon only, at its scale's positions; refused: an axis outside 0..2, a negative gap, a
  field not one per node.
- **CG6 (in part):** regions, zero gap and a field along x from a deck; the gap and the field change the image;
  refused: an unknown axis and a negative gap (usage), a 2D device.
- Looked at on the 3D MOSFET deck: the oxide and gate lifted off the silicon, and with n the inverted channel's
  surface (V_G = 5 V) between the pockets.
- **Mutation checks (CG7, Release): 9 of 9 caught** — boundary pieces put in region 0, interface pieces not turned,
  interface sides swapped, the default gap, the order reversed, the first region moved, the field ignored (caught
  once a field case was added), contacts left behind, `nitcad plot` ignoring `--field`.

**Commit 4 (filmstrips, transient currents, references, benchmarks; CG4, CG5, CG8), as built.** `scene.hpp`:
`FilmFrame`, `filmstrip_scenes` (one field view per frame on one scale fitted over every frame (CQ2), titled with
its label), `time_label` (three significant digits in s, ms, µs, ns, ps or fs). `image.hpp`: `tile_images` (rows of
`columns`, in order, on white). `nitcad plot --view filmstrip --field f [--every k] [--columns c] [--frames dir]`
on a transient result: at most 12 frames (every k-th snapshot from the first; by default every ⌈n/12⌉-th), each
drawn at the frame size (`--size`, else 640 × 480: a first version tiled the frames into `--size` and left each map
about 110 pixels wide), 4 to a row, each also written alone under `--frames`; a vector magnitude's ten-decade floor
taken over all frames. **Finding (fixed):** `nitcad plot` read a transient snapshot's fields and bands only; it now
keeps the snapshot's edge currents, rates and bias, so \|J\|, streamlines, arrows and R work on transient results.
**Arrows (amending commit 1):** `grid_arrows` and `surface_arrows` take an exclusion mask; a site whose nearest node
is masked gets no arrow, and `nitcad plot` masks the insulators for currents (the first 2D reference showed a row
of current arrows in the oxide, where the elements next to the interface interpolate the silicon's current).
- **CG4:** tiles in order, the empty place white, sizes checked; frames on one scale (1..5 from frames of 1..2 and
  3..5), titled, a given scale kept; time labels; a transient deck's filmstrip of n with every snapshot (frame count
  in the output, the strip's size, one PNG per frame), of \|J\| every 2nd in 2 columns, and current streamlines from
  the transient file; refused: a sweep's results, `--every 0`.
- **CG5:** WARP references (960 × 720 unless stated), each looked at before being kept: `mosfet2d_current_arrows`
  (\|J\| with 250 grid arrows, none in the oxide), `mosfet3d_efield_arrows` (E arrows on the exposed surface,
  colored by \|E\|), `mosfet3d_clip_n` (n with the quarter cut, the gate slab notched), `mosfet3d_exploded` (the oxide
  and its gate lifted), `mosfet2d_filmstrip_n` (n over the five bias points of the gate and drain sweep, 3 to a row of
  480 × 360 frames: the drain-side depletion grows with VD). 20 references in all pass.
- **CG8 (Release, other jobs running, so upper bounds):** the 3D MOSFET: 400 surface arrows 0.007 s, the quarter-cut
  surface 0.010 s, the exploded view 0.044 s; the 1,053,696-tetrahedron device: 400 arrows 0.61 s (the locator
  included), the quarter-cut surface 0.67 s (37,632 polygons), the exploded view 2.0 s.
- **Mutation checks (CG7, Release): 10 of 10 caught** — tiles transposed, the scale from the first frame only, the
  label left out of the title, the unit threshold, the decimals, the exclusion off, the snapshots' currents dropped,
  the default every-k, the default columns, the frames not written. (`nitcad plot`'s insulator mask for arrows is
  wiring the tests do not observe; the mask itself is tested, and the reference shows its effect.) Unit 25c in all:
  40 of 40.

**CG7 suites and numerics:** Release 22/22 and Debug 16/16 at each of the four commits (in clean worktrees), no
warnings under `/W4 /WX`, the Debug render suite's D3D12 debug layer clean; nothing outside `render`, `nitcad plot`,
tests, `CMakeLists.txt` and this file differs from `main` (`8256cf8`), so no numerical change. **Status:**
**accepted by the owner (2026-10-22, "accept")**, CQ1–CQ4 with it; ARCHITECTURE.md 6.14, 11 and 11.3 as built (CG9).
Push, PR and merge await the owner's authorization.
