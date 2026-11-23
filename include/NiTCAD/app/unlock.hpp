// Unlock shapes (Unit 26n, slice 2: DT5). No Win32: the dimensions panel's Unlock button calls it.
//
// A 3D deck the modeller cannot edit (a tensor or an extrude3d mesh) becomes a layered generate3d deck it can, by
// minimal edits through the deck text (every line outside [mesh], the regions' solid lines and the contacts' box
// lines kept byte for byte):
// - extrude3d: the modeller's own conversion (Unit 26l SQ5: read_model, convert_to_layered, write_model);
// - tensor (x, y and z): the domain is the mesh lines' ends as written (references kept), each region box stays a
//   box solid (references kept; a box reaching past the domain is cut to it), each contact's face and box become a
//   patch box on that face (references kept), max_size is the widest mesh-line spacing and interface_resolution
//   the largest the layout's rules allow (26m MT12's automatic s). The straddle lines go with the mesh lines.
// One-way: the app records the old text in the deck's undo. A refusal says why (an [adapt] section, which refines a
// tensor mesh; a layout the 3D generator refuses, in its words).
#pragma once

#include <expected>
#include <string>
#include <string_view>

namespace NiTCAD::app::unlock {

enum class Kind {
    none,      // not a 3D deck the modeller cannot edit (or it does not parse)
    tensor,    // a 3D tensor deck
    extruded,  // an extrude3d deck
};

[[nodiscard]] Kind kind_of(std::string_view deck);

struct Unlocked {
    std::string text;  // the layered deck
    std::string note;  // what changed, for the status bar
    double nodes = 0;  // the expected node count (MT12 (3)), 0 when not estimated
};

[[nodiscard]] std::expected<Unlocked, std::string> unlock(std::string_view deck);

}  // namespace NiTCAD::app::unlock
