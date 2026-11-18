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

}  // namespace NiTCAD::render
