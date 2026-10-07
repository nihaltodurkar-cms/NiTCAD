// Device description (ARCHITECTURE.md 6.3): the problem a solve runs on, as plain data.
//
// A device is a mesh, a lattice temperature, regions (each a name and a material), the region of
// every node, the donor and acceptor concentrations of every node, and the contacts. It holds no
// solver state and is built without any solver present. It is validated once, at construction, so
// the layers above can rely on it.
//
// Doping is given as separate donor and acceptor concentrations because the equations need two
// different combinations: the net doping N_D - N_A (Poisson charge) and the total ionised impurity
// N_D + N_A (mobility and lifetime). The legacy took the net doping and, by default, |net| as the
// total; that is the case where each node carries one dopant type.
//
// Several regions with different materials are allowed by the data model; the assemblers treat an
// edge between two materials as a heterointerface (Unit 15, assemble/scaled_device.hpp), across
// which carriers drift and diffuse. A declared Interface between two regions may choose another
// transport law for the edges joining them (ARCHITECTURE.md principle 9: interfaces carry their own
// models). A graded composition, described as many regions, needs no declarations: its steps are
// drift-diffusion.
#pragma once

#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/contact.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::device {

using RegionId = std::int32_t;

struct Region {
    std::string name;
    physics::Semiconductor material;
};

// The carrier transport across the edges joining two regions.
enum class InterfaceTransport : std::uint8_t {
    drift_diffusion,      // Scharfetter-Gummel with the band offsets (the default for any edge)
    thermionic_emission,  // emission-limited fluxes (assemble/thermionic_flux.hpp)
};

// A declared interface between two regions (order does not matter).
struct Interface {
    std::string region_a;
    std::string region_b;
    InterfaceTransport transport = InterfaceTransport::drift_diffusion;
};

struct DeviceDescription {
    mesh::Mesh mesh;
    double temperature_K;               // lattice temperature, uniform (no self-heating)
    std::vector<Region> regions;
    std::vector<RegionId> node_region;  // one per mesh node, indexing regions
    std::vector<double> donors;         // N_D per mesh node [cm^-3]
    std::vector<double> acceptors;      // N_A per mesh node [cm^-3]
    std::vector<Contact> contacts;
    std::vector<Interface> interfaces;  // declared interfaces; may be empty
};

class Device {
public:
    // Validates the description. Every error is invalid_input; the context index names the
    // offending region, node or contact, as the message says:
    // - no regions; a region name empty or repeated; the temperature rejected by
    //   physics::check_temperature for a region's material (the message names the region);
    // - node_region, donors or acceptors not one per mesh node; a node's region out of range; a
    //   region with no nodes;
    // - a donor or acceptor concentration that is not finite and >= 0 (the value is in context);
    // - a contact name empty or repeated, or of unknown kind; a contact with no nodes, nodes out
    //   of range or not strictly increasing, a node on no mesh boundary patch, or a node shared
    //   with another contact;
    // - a gate contact whose stack is invalid (contact.hpp: oxide thickness or permittivity not
    //   finite and positive, an unknown electrode, a metal work function not finite and positive,
    //   a fixed charge not finite), whose boundary patch does not exist, or with a node not on
    //   that patch;
    // - a floating region: a connected part of the mesh graph with no ohmic contact node. Its
    //   potential or carrier densities would be undetermined (the Poisson or continuity matrices
    //   singular), the case the linear solver's pivot-ratio heuristic otherwise has to catch
    //   (6.10); a gate does not anchor the carriers. The index is the lowest node of that part;
    // - an interface naming an unknown region, the same region twice, a pair already declared (in
    //   either order), a pair no mesh edge joins, or an unknown transport; the index is the
    //   interface.
    [[nodiscard]] static std::expected<Device, base::Error> create(DeviceDescription description);

    [[nodiscard]] const mesh::Mesh& mesh() const noexcept { return d_.mesh; }
    [[nodiscard]] double temperature_K() const noexcept { return d_.temperature_K; }
    [[nodiscard]] std::span<const Region> regions() const noexcept { return d_.regions; }
    [[nodiscard]] std::span<const RegionId> node_region() const noexcept { return d_.node_region; }
    [[nodiscard]] std::span<const double> donors() const noexcept { return d_.donors; }
    [[nodiscard]] std::span<const double> acceptors() const noexcept { return d_.acceptors; }
    [[nodiscard]] std::span<const Contact> contacts() const noexcept { return d_.contacts; }
    [[nodiscard]] std::span<const Interface> interfaces() const noexcept { return d_.interfaces; }
    // The transport across an edge joining regions a and b: the declared interface's, else
    // drift-diffusion. Precondition (NITCAD_EXPECTS): both are valid region ids.
    [[nodiscard]] InterfaceTransport transport(RegionId a, RegionId b) const;
    // The contact with this name, or nullptr.
    [[nodiscard]] const Contact* find_contact(std::string_view name) const noexcept;

    // Per-node conveniences. Precondition (NITCAD_EXPECTS): node is a valid node id.
    [[nodiscard]] double net_doping(mesh::NodeId node) const;      // N_D - N_A [cm^-3]
    [[nodiscard]] double total_impurity(mesh::NodeId node) const;  // N_D + N_A [cm^-3]
    [[nodiscard]] const physics::Semiconductor& material(mesh::NodeId node) const;

private:
    explicit Device(DeviceDescription description) noexcept : d_(std::move(description)) {}

    DeviceDescription d_;
};

}  // namespace NiTCAD::device
