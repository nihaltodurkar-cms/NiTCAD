<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Units 25 and 26 — visualization and the Windows application (proposed and authorized 2026-10-17; on `render/visualization`)

**Owner request (2026-10-17):** a scope for the "Visualization, GUI, batch" row of ARCHITECTURE.md 11.3 (Units 25
and 26, both "target"), presented as an HTML page; then the target look given by the owner's reference grid of
twelve 3D MOSFET views ("TCAD Semiconductor Visualisation Grid"). The proposal, with mockups of the window and of
the twelve views, is the private page https://claude.ai/artifact/A3K7DSo1FZMjsm1aew5htX (version 3); this entry is
its record.

**Why now (audited on `main` at `913dac0`).**
- Every result type has a file since Unit 42 and every analysis runs from a deck since Unit 43; the only output a
  user sees is `nitcad`'s text and CSV tables.
- The roadmap reserves `render` (L7: base, mesh, results, analysis) and `app` (L8: everything except physics and
  assemble). 6.9 fixes the app's threading: every solve on a worker thread, progress posted to the UI thread.
- On this workstation: Windows SDK 10.0.26100 with `d3d12.h`, `d2d1_3.h`, `d3d11on12.h`, `dxgi1_6.h`,
  `wincodec.h` and `dxc.exe`; WARP (`d3d10warp.dll`); an NVIDIA GeForce RTX 5060 Ti. No new vcpkg dependency.
- A mesh is a node and edge graph without cells: surfaces come from a generated mesh's simplices or a tensor
  device's `TensorCells`. A sweep file holds fields, not the mesh: the viewer pairs it with the deck's mesh file.

**Design: Unit 25 (`render`).**

| # | Choice |
|---|---|
| VR1 | A CPU scene (vertices, values, colors, line strips, labels) built from a mesh and a result, separate from drawing; drawing never reads results. |
| VR2 | Per dimension: 1D line plots; 2D field maps on triangles (simplices, or tensor cells split in two) with edges, region outlines and contacts; 3D the reference grid's views: region-colored solid view with contacts and electrodes, any field on the boundary surface, the tetrahedral mesh with regions, slice planes (XY, XZ, YZ), isosurfaces, current streamlines, transparency. 3D tensor cells are split into six tetrahedra. A mesh without simplices or cells is drawn as nodes and edges only. |
| VR3 | Fields: potential, n, p, |E|, |J|, band edges, quasi-Fermi levels, temperature, doping; linear and log scales; colorbars with units; contours by marching triangles; isosurfaces by marching tetrahedra; streamlines of J traced through tetrahedra (RK4, seeded on a contact or a line). |
| VR4 | Plots (I-V, C-V, transient waveforms, admittance against frequency, band diagrams) from `analysis` curves, with nice ticks, log axes, units and Unit 24's resolution and extrapolation flags. |
| VR5 | Drawing: triangles and lines on Direct3D 12; text, axes and legends on Direct2D through D3D11On12; off-screen targets read back and saved as PNG through WIC; a hardware device when present, else WARP. |
| VR6 | `nitcad plot <deck or result file> ...`: fields, views and curves to PNG; a result file is paired with its mesh (deck outputs or `--mesh`), node counts checked. |
| VR7 | Out of scope: windows (Unit 26), animation, ray-marched volume rendering, VTK export, printing. |
| VR8 | The reference's look: lit, shaded solids; fixed region colors (substrate beige, n⁺ red, p channel blue, oxide pale translucent, metal dark grey); an electrode drawn as a slab of nominal thickness on its contact patch (a gate is a contact in NiTCAD, not a meshed region); weighted blended order-independent transparency; an axis triad in every view. |
| VR9 | Viewing reconstructions in `render`: |E| per tetrahedron from the potential's gradient, J per tetrahedron from the edge currents (least squares over its six edges); nothing added to the solver or to files. Recombination is not in result files and `render` may not include `physics`, so the reference's R view waits for a unit that stores the rate in results. |

**Design: Unit 26 (`app`).** VA1 one main window (Win32, per-monitor DPI v2: deck pane, views, run bar, log); VA2
the deck pane (Rich Edit, line numbers, Check puts `deck::build`'s errors at their lines); VA3 runs on a worker
thread, progress by posted message, Stop through the stop token, completed parts kept; VA4 Unit 25's scenes in a
swap chain with pan, zoom, rotate and node picking; VA5 a results browser (run record, provenance, convergence,
Unit 24's figures); VA6 logic in plain C++ under test, a thin window layer, `--self-test` for CI on WARP; VA7 out of
scope: a graphical device editor, projects, remote runs, plugins, localization, an installer.

**Owner rulings (2026-10-17):**

| # | Ruling |
|---|---|
| VQ1 | D3D11On12 for Direct2D on Direct3D 12. |
| VQ2 | HLSL compiled at build time with the SDK's `dxc`, bytecode embedded. |
| VQ3 | Units in sequence: 25a, then 25b, then 26. |
| VQ4 | Exact tests of the CPU scene, plus WARP image regression tests against committed references. |
| VQ5 | Turbo as the default colormap; viridis, cividis and a diverging map also supported. |
| VQ6 | Unit 44 (iterative solver) first if it is approved; then Unit 25, then Unit 26. |
| VQ7 | Unit 25 split: 25a, the core rendering (scene, surfaces, solid and mesh views, slices, transparency, the image pipeline) and `nitcad plot`; 25b, derived |E| and |J|, isosurfaces and streamlines. |

**Owner conditions (2026-10-17):** numerical behaviour and the architecture's boundaries are preserved; Release
and Debug tests and warning-free `/W4 /WX` builds; scientific validation and performance benchmarks; only the
currently approved unit is implemented, and each is accepted before the next starts; nothing is pushed or merged
without authorization.

**Acceptance gates.**

| # | Gate | Unit |
|---|---|---|
| VG1 | Layering: `render` includes only base, mesh, results and analysis (a header-boundary test); the assembly traces are byte-identical to `main`'s (no numerical change). | 25a |
| VG2 | Geometry: triangles of tensor and generated 2D meshes cover the domain exactly once; 3D boundary surfaces are closed; tetrahedra from tensor cells fill each cell; slice polygons match analytic sections. | 25a |
| VG3 | Values: colormap tables and endpoints; log scales exact at decades; contours of a linear field at the exact positions; ticks of 1, 2 or 5 times a power of ten. | 25a |
| VG4 | Images: WARP renders of a band diagram, a 2D potential map, an I-V plot and the 3D MOSFET's solid, surface (n, p, φ), mesh, slice and transparent views match committed references within a stated per-pixel tolerance; in 25b also |E|, |J| with streamlines and isosurfaces. The D3D12 debug layer reports no errors in Debug. | 25a, 25b |
| VG5 | `nitcad plot` writes each image from decks and from result files and refuses a mismatched mesh (exit 1). | 25a |
| VG6 | Release and Debug suites at every commit, no warnings under `/W4 /WX`; hand mutants in the scene algorithms are caught (run in Release). | 25a, 25b |
| VG7 | Benchmarks: scene building and drawing of 10⁵ and 10⁶ triangles on WARP and on the GPU, recorded; the GPU budget is 60 frames per second at 10⁶ triangles. | 25a |
| VG8 | Scientific validation: marching tetrahedra on a linear field gives the exact plane, closed and oriented; on a sphere's distance field the enclosed volume converges at second order; streamlines of a uniform field are straight and of a rotational field close within the integrator's tolerance; |E| of a linear potential and J of a uniform current are exact per tetrahedron. | 25b |
| WG1 | No solve on the UI thread: the window keeps pumping messages (within 100 ms) during a long run; Stop keeps the completed part. | 26 |
| WG2 | A deck run in the app is bit-identical to `nitcad run` on it, run identity included. | 26 |
| WG3 | The app's logic (run state machine, documents, view settings) under unit tests; deck errors land on their lines. | 26 |
| WG4 | `--self-test` passes on WARP and on the GPU and writes the same image as `nitcad plot` for the same view. | 26 |
| WG5 | Per-monitor DPI re-layout, keyboard access to every command, a clean D3D12 debug layer, no leaked GPU objects. | 26 |
| WG6 | Suites as VG6; nothing outside VA1–VA6; ARCHITECTURE.md updated to "as built" after acceptance. | 26 |

**Commits (planned):** 25a·1 the CPU scene (VG2, VG3); 25a·2 off-screen D3D12 and Direct2D drawing, lighting, region
colors, electrode slabs, transparency, PNG, references and the benchmark (VG4, VG7); 25a·3 `nitcad plot` (VG5);
25b derived fields, isosurfaces, streamlines (VG8, VG4); 26·1 the app's logic (WG1, WG3); 26·2 the window (WG2,
WG4, WG5). Estimated sizes: 25a about 3,900 lines, 25b about 1,700, 26 about 2,900, tests included.

**Risks:** GPU images differ between vendors and drivers (image gates on WARP only); D3D11On12's per-frame cost;
large 3D results (single-precision GPU copies, scales from the doubles); `dxc` and WARP on the CI runner are
unverified (checked in 25a's first commits); screen-reader support for custom views is out of scope.

**Detailed plan (page versions 4–6, 2026-10-17).** At the owner's request the page planned the application in
more detail and took the look of the owner's workspace reference (a light 3D workspace: menu bar, tool bar, scene
tree, viewport with overlay cards, metrics strip, inspector, status footer). VA1 now describes that workspace; VA8
adds the viewport tools and overlays (select, orbit, pan, zoom, fit; axonometric and perspective; view, field,
colorbar, probe and orientation overlays; probe and measure in 26, slice in 25a, clip in 25c). Figures 3–8 plan
the start page, the deck editor, the run monitor, the results workspace, cut lines with band diagrams and run
comparison; sixteen more views (13–28) are assigned to units. Sizes: 25a about 4,100 lines, 25b 1,700, 25c 1,600,
26 3,300, 26b 2,600, 45 500.

**Rulings VQ8–VQ14 (2026-10-17).** The owner delegated these to the recommendation best for a commercial TCAD
product ("for decisions choose whats best for Commercial TCAD"); each is the page's recommendation.

| # | Ruling | Why, for a commercial product |
|---|---|---|
| VQ8 | `render` may depend on `device`, read-only: regions, material names and the insulator flag, contacts, doping. No solver, assembler or physics model. VG1 becomes: `render` includes only base, mesh, results, analysis and device. ARCHITECTURE.md 3 changes with 25a. | Region colors, electrodes and doping are what users look at first; a copied "view description" in results would be a second source of truth for the same geometry. |
| VQ9 | Unit 25c (vector glyphs, clipping box, exploded view, transient filmstrips) after 25b, each view droppable. | Vectors and clipping are expected in every commercial viewer; independent views keep the risk small. |
| VQ10 | Unit 26b (live deck checking and completion, structure preview, cut lines with band diagrams, run comparison and difference maps, mesh quality) after 26. | These make the application an analysis tool rather than a viewer; 26 keeps its size. |
| VQ11 | A top-level `examples/` directory with four decks (1D diode, 2D MOSFET, 3D MOSFET, MOS-C), each run by a test. | A first-run user needs working examples; tested examples never rot. |
| VQ12 | Unit 45 (recombination and generation rates in results and files; the mesh digest in sweep, transient and small-signal files) before 25b. | Unlocks the R view and exact result–mesh pairing instead of a node-count check. |
| VQ13 | The Device, Mesh, Physics and Simulation menus open the deck at that section and list its keys for insertion; the deck stays the only source of truth. Form editors that write the deck may be a later unit. | Reproducible, diffable, scriptable input is what production flows depend on; a form editor that writes the deck can come later without a second model. |
| VQ14 | Turbo stays the default (VQ5); the reference's softer five-color map is added as a named option, `soft`. | Turbo is the familiar high-contrast default; `soft` matches the reference where wanted. |

**Status:** authorized with VQ1–VQ14. Unit 44 was implemented first (VQ6) and merged (`b518c99`, PR #54); Unit 25a
started on the owner's "go 25a" (2026-10-21), on `render/visualization` rebased onto `b518c99`; commits 25a·1–3 done;
**accepted by the owner (2026-10-21, "accepted")**; ARCHITECTURE.md 6.14, 3, 11 and 11.3 as built. Then 45, 25b, 25c,
26, 26b, each approved after the one before it is accepted.

**Commit 25a·1 (the CPU scene; VG1–VG3, VG6), as built.** The `render` layer (`include/NiTCAD/render/`,
`src/render/`; CMake target `NiTCAD::render`, linked to base, mesh, results, analysis and device):
- `color.hpp`: Turbo, viridis, cividis (matplotlib 3.9.2's 256-entry tables, extracted with awk), Moreland's
  cool-warm (33 points) and `soft` (VQ14: the reference names deep blue, teal, green, gold and muted red without
  values; chosen as #253f8f, #2a9d9c, #6dbb5a, #e3b23c, #c0504d), interpolated linearly in sRGB; `ColorScale` linear,
  log10 and signed log10 (sign(v) log10(1 + |v| / floor), for doping); ticks of 1, 2 or 5 times a power of ten,
  decades on log axes (minor 2..9 on axes of at most three decades), labels with U+2212 and superscript powers.
- `geometry.hpp`: elements from `TensorCells` (2 triangles; 6 Kuhn tetrahedra, conforming between cells) or a
  `SimplexMesh` (reoriented counterclockwise), else the graph; 3D boundary faces (outward); slices (a facet lying
  in the plane taken once, from its positive side); **regions drawn by control volume**: in the box method the
  interface between two regions runs between their nodes (no generated node lies on one), so a tensor cell is
  split at its mid-planes (exactly the box-method control volumes) and a simplex by its planes of equal
  barycentric weight; pieces, interfaces and both rules by convex-polygon clipping; marching-triangle contours (a
  vertex at the level counts as above it).
- `scene.hpp`: `FieldScene` (triangles with per-vertex colormap positions or flat colors, line segments, colorbar,
  camera, bounds, in µm) and `PlotScene` (axes, series, flags); views `solid_view`, `mesh_view`, `field_view`
  (2D maps; 3D on the boundary surface, with the semiconductor seen through translucent insulators at their
  interfaces), `slice_view`; contact slabs of nominal thickness (3% of the device: a lumped gate as oxide under
  metal, an electrode as metal, an ohmic contact as a thin pad); the VR8 palette by doping class (n⁺ ≥ 1e19 red,
  n salmon, p channel 1e17–1e19 blue, p⁺ dark blue, substrate beige, insulators pale at alpha 0.35, metal dark
  grey); `field_values` for ψ, n, p, the bands and quasi-Fermi levels, T and doping; `curve_plot` (Unit 24 curves;
  a log axis of all-negative values plots |y|; points at or below their resolution flagged), `profile_plot`,
  `band_plot` (1D; insulators as gaps).
- **Gates:** `tests/render` (26 cases): VG2 — 2D tensor and generated triangles cover the domain exactly once
  (no directed edge twice); 3D tensor and generated boundaries closed and outward (enclosed volume and area to
  1e-12); Kuhn tetrahedra fill each cell; slices match the analytic sections (x = 0.73, a grid plane, the hexagon
  x + y + z = 1.5) with a linear field exact on them; the control-volume pieces of a tensor mesh equal its control
  volumes to 1e-14; tensor interfaces lie halfway between the nodes; VG3 — table endpoints, log positions exact at
  decades, ticks and labels; scenes — region areas per material equal the control volumes, field positions exact,
  contours on their levels, winding consistent with every normal, the silicon under an oxide facing the oxide
  above and below. VG1: ctest `render_layering` (a CMake script, no Python) checks that `render` includes only
  base, mesh, results, analysis and device; nothing outside `render`, its tests and CMakeLists.txt changed.
- **Mutation checks (VG6, Release): 18 of 18 caught** (Kuhn, simplex and face orientation, the slice facet rule,
  the tensor and barycentric owner rules, the interface position, the contour tie rule, the table step, the ticks'
  2, the scale's clamp, the fan's orientation, the interface normal, insulator transparency, field and slice
  positions, the doping threshold, a forbidden include). A degenerate-triangle check was found dead (the vertex
  merge before it removes every case) and removed. Release 22/22 and Debug 16/16 (with the new ctests `render` and
  `render_layering`), no warnings under `/W4 /WX`.

**Commit 25a·2 (drawing; VG4, VG6, VG7), as built.** `image.hpp`: a `Renderer` on a Direct3D 12 device (the first
hardware adapter, else WARP) draws a `FieldScene` or a `PlotScene` off screen into an `Image`; `save_png`,
`load_png` (WIC) and `compare`.
- Triangles and lines on Direct3D 12 at twice the image's size, averaged 2 × 2; field triangles look the colormap up
  per pixel from the interpolated position (a 256-entry table resampled from `sample`); one directional light
  (from the upper left over the viewer's shoulder) and ambient 0.38, two-sided; translucent triangles by weighted
  blended OIT (McGuire and Bavoil); lines as screen-space quads of their width, depth-tested or on top. Text,
  colorbars, 2D axes in µm, the axis triad and plots on Direct2D and DirectWrite (Segoe UI, grayscale text
  antialiasing) drawn into the same target through D3D11On12 (VQ1). HLSL compiled at build time by the SDK's `dxc`
  into headers of bytecode (VQ2; CMake finds `dxc` and fails without it).
- **The camera convention (changes 25a·1's):** +y up, the device's vertical axis, as in 2D; azimuth about +y from +z
  toward +x, elevation above the xz plane; 2D looks down −z. The default 3D view is azimuth 35°, elevation 25°.
  Drawing fits the scene's bounds to the image's aspect (92% of the view).
- **Debug layer:** the one message the pipeline raises, `REFLECTSHAREDPROPERTIES_INVALIDOBJECT` from D3D11On12's
  `CreateWrappedResource` probing a resource Direct3D 12 created, is filtered by its exact ID; every other error
  counts, and the tests require none (Debug builds run the image tests with the layer).
- **Errors:** a device failure is `resource_exhausted` with the call and its HRESULT (there is no device error code in
  `base`; a new one was not added in this unit).
- **VG4 references** (`tests/render/references/`, 11 images, 960 × 720, WARP): a 2D MOSFET's potential map with
  contours and its regions with the mesh, a band diagram, an I-V plot with flagged points, and a 3D MOSFET's solid,
  n, p and ψ surfaces, mesh, two slices, and an opaque perspective view; the MOSFETs are solved in the test (21 × 21
  and 21 × 21 × 7 nodes, Vg = 1 V, Vd = 0.5 V). Each was looked at before it was committed; the first renders
  showed the views mirrored (the camera's right vector reversed), fixed before any reference was kept. A render
  matches when at most 0.5% of its pixels differ by more than 24 in a channel; `NITCAD_UPDATE_REFERENCES=1` rewrites
  them. Also exact checks: a flat triangle's color, a field's colormap color, image sizes, the PNG round trip, and
  order-independent transparency giving the same image in either order.
- **VG7 benchmark** (`[.benchmark]`, 1280 × 960, the boundary of a thin tensor box):

  | Triangles | Adapter | Image (draw: upload, passes, overlays, readback) | Redraw on resident geometry |
  |---|---|---|---|
  | 6.8e4 | RTX 5060 Ti | 13.6 ms (74/s) | 0.20 ms |
  | 1.0e6 | RTX 5060 Ti | 146 ms (6.8/s) | 0.81 ms (1236 frames/s) |
  | 6.8e4 | WARP | 29.0 ms | 13.8 ms (72/s) |
  | 1.0e6 | WARP | 302 ms | 88.1 ms (11/s) |

  The 60 frames per second budget at 10⁶ triangles is met on the GPU by the passes a window repeats (0.81 ms); a
  whole image costs 146 ms, nearly all of it the per-image upload of 132 MB of vertices and the readback, which a
  window does not repeat. Scene building for 10⁶ triangles takes 1.3 s on the CPU.
- **Mutation checks (VG6, Release): 10 of 10 caught** — the camera's handedness, the fit margin, the readback and
  PNG channel order, the field flag, a line quad's corner, the depth test, and in the shaders the lighting, the OIT
  revealage and the 2 × 2 average.
- ARCHITECTURE.md 3: `render`'s row of allowed dependencies adds `device`, read-only (VQ8).
- **Suites:** Release 22/22 and Debug 16/16, no warnings; `nitcad_render_test` 30 test cases (+4). The Debug `render`
  ctest takes 641 s (the two MOSFETs solved unoptimized, then drawn on WARP with the debug layer), Release 12 s.

**Commit 25a·3 (`nitcad plot`; VG5), as built.** `nitcad plot <deck | result file> [options] --out <image.png>`
(`src/cli/plot.cpp`; `nitcad` now links `render`):
- **Input:** a deck (its device and mesh rebuilt, the result file its analysis writes read: `sweep`, `transient`,
  `small_signal`, or the equilibrium's `checkpoint`), or a `.nitsweep`, `.nittran` or `.nitac` file with `--device
  <file.nitdev>` and, for a generated mesh, `--mesh <file.nitmesh>` (its simplices; refused unless its graph digest is
  the device mesh's). A result whose fields are not one value per node of the device's mesh is refused (exit 1):
  the node-count pairing of VR6, until Unit 45 stores the mesh digest in result files.
- **Views:** `solid`, `mesh`, `field` (a 1D device's is a profile plot), `slice` (`--slice x|y|z=<um>`, repeatable),
  `profile`, `bands` (1D), and the curves `iv` (a sweep, `analysis::current_curve`), `cv` (a small-signal run,
  `capacitance_curve`) and `current` (a transient's terminal currents against time). Options: `--field` (the keys of
  `field_info`), `--point`, `--scale`, `--range`, `--colormap`, `--contours`, `--edges`, `--opaque`, `--azimuth`,
  `--elevation`, `--perspective`, `--contact`, `--sense`, `--frequency`, `--log`, `--size` (1280 × 960), `--adapter`.
  Regions and doping need no results; other views before the run fail with "(run the deck first)".
- **Exit status:** 0 written, 1 failure (the message on stderr), 64 usage (the option list on stderr).
- **Gates** (`tests/deck/run_test.cpp`, the executable started as the other CLI tests do, `--adapter warp`): images
  of the expected size from a deck (regions before the run; n with contours, the mesh, doping and I-V after it),
  from a result file with its device (a chosen point), and from a generated mesh's result with `--mesh`; refusals
  with exit 1 of a views-before-run request, a point out of range, another deck's device (node counts named) and
  another deck's mesh file; usage errors with exit 64. Checked by hand on a 2D pn diode: n across the junction and an
  I-V of about 60 mV per decade. (The first try drew a uniform n: the demo deck's profiles had their depth axis
  along x, so their x windows were ignored and the doping compensated to zero, which the plot showed as it was.)
- **Mutation checks (VG6, Release): 7 of 7 caught** — the node-count and mesh-digest refusals, the point's range
  (first caught only by the status, through undefined behaviour; the test now checks the message), the results a view
  needs, `--out` and `--device` required, the `--point` option.
- **VG1, assembly traces:** nothing in Unit 25a changes `assemble`, `solve`, `deck` or `io` (render, its tests, the
  CLI's plot command and CMakeLists.txt only), so every existing suite runs the code `main` runs; the suites pass
  unchanged. Release 22/22 and Debug 16/16, no warnings.
