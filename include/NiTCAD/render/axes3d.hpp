// Scale and orientation overlays of the enhanced 3D look (Unit 26h; DECISIONS.md JD1-JD5, JQ1-JQ6).
//
// axis_box lays out the axis box of a 3D FieldScene drawn at a size: the edges of the three faces of the
// scene's bounds turned away from the camera (ParaView's "grid axes", JQ2 (a)), and along one edge per axis
// 1-2-5 ticks in um with their labels. Each axis takes one of the box's outline edges parallel to it, both
// on the far side from the box's centre: x and z the one nearer the viewer, y the one further left on the
// screen (JQ1 (a); JD1's own example, "the front y edge on the left"). A face seen exactly edge-on counts as
// turned away on its lo side. An axis whose extent is zero, or whose edge is under 24 px on the screen
// (looking along it), has no ticks. Labels nearer a kept label (of any axis) than their own extent plus 6 px
// are dropped (widths estimated at 6.5 px per character of 11.5 px Segoe UI). The renderer draws the edges depth-tested
// and the ticks and labels over the image; nothing here touches the GPU.
#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::render {

struct AxisLabel {
    Float3 position{};  // the tick on its edge, um
    float x = 0.0f;     // its pixel
    float y = 0.0f;
    float label_x = 0.0f;  // where the label is anchored (centred), pixels
    float label_y = 0.0f;
    double value = 0.0;  // um
    std::string text;
    bool kept = true;  // false: dropped as overlapping a kept label
};

struct BoxAxis {
    int axis = 0;                 // 0 x, 1 y, 2 z
    std::array<Float3, 2> edge{};  // the edge it lies on, um
    std::array<float, 2> outward{};  // unit, on the screen, away from the box's centre
    std::vector<AxisLabel> ticks;
    std::string title;  // "x (µm)"
    float title_x = 0.0f;
    float title_y = 0.0f;
};

struct AxisBox {
    // For each axis, whether the face turned away from the camera is the hi side (else the lo side).
    std::array<bool, 3> away_hi{};
    std::vector<std::array<Float3, 2>> edges;  // the away faces' edges, each once (9 for a solid box)
    std::vector<BoxAxis> axes;                 // the axes that have ticks
};

// Errors none: a 2D scene gives an empty box.
[[nodiscard]] AxisBox axis_box(const FieldScene& scene, int width, int height);

// --- the view cube (Unit 26h, JD2, JD3) -------------------------------------------------------------
//
// view_cube lays out the enhanced look's orientation cube, which replaces the axis triad in the view's
// lower-left corner: a cube 72 px across its diagonal, turned with the camera (orthographic, whatever the
// scene's projection), its faces toward the eye with their names (Front +z, Back -z, Right +x, Left -x,
// Top +y, Bottom -y), shaded by the key light, and x, y and z arrows from its corner nearest the lower left
// of the screen, along the triad's directions. view_cube_hit gives the camera a click there asks for: a
// corner of a face shown (within 7 px) its isometric view, else a face shown its named view, else nothing.
// Both are empty for a classic, 2D or triad-less scene.

struct CameraPreset {
    std::string name;
    float azimuth_deg = 0.0f;
    float elevation_deg = 0.0f;
};

// The named views (JD4): Front, Back, Left, Right, Top, Bottom, and Isometric (the default camera).
[[nodiscard]] const std::vector<CameraPreset>& named_views();

struct CubeFace {
    std::string label;  // "Front", ...
    int axis = 0;       // its normal: +-x, y or z
    int sign = 1;
    std::array<std::array<float, 2>, 4> polygon{};  // pixels, in order around the face
    float shade = 1.0f;                             // 0.62 to 1, by the key light
    CameraPreset preset;
};

struct CubeCorner {
    std::array<float, 2> at{};  // pixels
    std::array<int, 3> signs{};  // which corner: -1 or +1 along x, y, z
    CameraPreset preset;         // the isometric view from it
};

struct CubeArrow {
    std::string name;  // "x", "y", "z"
    std::array<float, 2> from{};
    std::array<float, 2> to{};  // empty (from == to) when the axis points at the eye
};

struct ViewCube {
    std::array<float, 2> centre{};
    float edge_px = 0.0f;
    std::vector<CubeFace> faces;      // the faces toward the eye
    std::vector<CubeCorner> corners;  // the corners of those faces
    std::vector<CubeArrow> arrows;    // x, y, z
};

[[nodiscard]] ViewCube view_cube(const FieldScene& scene, int width, int height);
[[nodiscard]] std::optional<CameraPreset> view_cube_hit(const FieldScene& scene, int width, int height, float x,
                                                        float y);

}  // namespace NiTCAD::render
