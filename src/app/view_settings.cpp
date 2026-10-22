// View settings into scenes (view_settings.hpp). Each view is built with the calls and choices of
// `nitcad plot` (src/cli/plot.cpp) for the same options, so the window and the command draw the same
// picture (WG4); the app does not include the command's code (AQ4), and the WG4 test compares the two
// images.
#include "NiTCAD/app/view_settings.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <format>
#include <utility>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/app/compare.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/render/compare.hpp"
#include "NiTCAD/render/quality.hpp"
#include "NiTCAD/render/streamlines.hpp"
#include "NiTCAD/render/vectors.hpp"

namespace NiTCAD::app {

namespace {

using namespace NiTCAD::render;

constexpr std::array<std::pair<ViewKind, std::string_view>, 14> view_names{{{ViewKind::solid, "solid"},
                                                                          {ViewKind::mesh, "mesh"},
                                                                          {ViewKind::field, "field"},
                                                                          {ViewKind::slice, "slice"},
                                                                          {ViewKind::iso, "iso"},
                                                                          {ViewKind::bands, "bands"},
                                                                          {ViewKind::profile, "profile"},
                                                                          {ViewKind::iv, "iv"},
                                                                          {ViewKind::cv, "cv"},
                                                                          {ViewKind::current, "current"},
                                                                          {ViewKind::quality, "quality"},
                                                                          {ViewKind::quality_histogram, "quality-histogram"},
                                                                          {ViewKind::difference, "difference"},
                                                                          {ViewKind::exploded, "exploded"}}};

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

// A number as `nitcad plot` reads its options: the whole text, finite.
std::optional<double> number(std::string_view s) {
    double v = 0.0;
    const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size() || !std::isfinite(v)) return std::nullopt;
    return v;
}

// A segment "<x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>]" in um, as its ends in cm, each end with exactly `dimension`
// coordinates (as `nitcad plot --seeds` reads it).
std::optional<std::pair<Point, Point>> segment_um(std::string_view s, int dimension) {
    const auto point_of = [&](std::string_view t) -> std::optional<Point> {
        Point q{};
        std::size_t k = 0;
        for (std::string_view rest = t;; ++k) {
            const auto comma = rest.find(',');
            const auto x = number(rest.substr(0, comma));
            if (!x || k >= 3) return std::nullopt;
            q[k] = *x * 1e-4;
            if (comma == std::string_view::npos) break;
            rest = rest.substr(comma + 1);
        }
        if (k + 1 != static_cast<std::size_t>(dimension)) return std::nullopt;
        return q;
    };
    const auto colon = s.find(':');
    if (colon == std::string_view::npos) return std::nullopt;
    const auto a = point_of(s.substr(0, colon));
    const auto b = point_of(s.substr(colon + 1));
    if (!a || !b) return std::nullopt;
    return std::pair{*a, *b};
}

std::optional<VectorField> vector_field(const std::string& name) {
    if (name == "E") return VectorField::electric;
    if (name == "J") return VectorField::current;
    if (name == "Jn") return VectorField::electron_current;
    if (name == "Jp") return VectorField::hole_current;
    return std::nullopt;
}

std::string current_unit(int dimension) { return dimension == 1 ? "A/cm\xC2\xB2" : dimension == 2 ? "A/cm" : "A"; }

std::size_t contact_index(const device::Device& d, const std::string& name, std::size_t fallback) {
    for (std::size_t k = 0; k < d.contacts().size(); ++k) {
        if (d.contacts()[k].name == name) return k;
    }
    return fallback;
}

}  // namespace

std::string_view to_string(ViewKind v) noexcept {
    for (const auto& [k, name] : view_names) {
        if (k == v) return name;
    }
    return "?";
}

std::optional<ViewKind> find_view(std::string_view name) noexcept {
    for (const auto& [k, n] : view_names) {
        if (n == name) return k;
    }
    return std::nullopt;
}

bool is_plot(ViewKind v) noexcept {
    return v == ViewKind::bands || v == ViewKind::profile || v == ViewKind::iv || v == ViewKind::cv ||
           v == ViewKind::current || v == ViewKind::quality_histogram;
}

bool takes_vectors(ViewKind v) noexcept { return v == ViewKind::field || v == ViewKind::slice || v == ViewKind::iso; }

std::expected<Scene, base::Error> build_scene(const Workspace& w, const ViewSettings& settings, const Workspace* b) {
    // Arrows and streamlines on the field views only (Unit 26c, YD1).
    const bool drop = !takes_vectors(settings.kind) && (!settings.vectors.empty() || !settings.streamlines.empty());
    ViewSettings kept;
    if (drop) {
        kept = settings;
        kept.vectors.clear();
        kept.streamlines.clear();
    }
    const ViewSettings& o = drop ? kept : settings;
    if (!w.device || !w.geometry) return std::unexpected(invalid("nothing is loaded"));
    const device::Device& device = *w.device;
    const Geometry& geometry = *w.geometry;
    const mesh::Mesh& m = device.mesh();
    const auto& points = w.points;

    const auto pick = [&]() -> std::expected<const results::BiasPoint*, base::Error> {
        if (points.empty()) return std::unexpected(invalid("the results hold no field to plot"));
        const std::size_t k = o.point.value_or(points.size() - 1);
        if (k >= points.size()) {
            return std::unexpected(invalid("--point " + std::to_string(k) + " is out of range (" +
                                           std::to_string(points.size()) + " points)"));
        }
        const auto& p = points[k];
        if (p.fields.potential_V.size() != m.node_count()) {
            return std::unexpected(invalid("the result has " + std::to_string(p.fields.potential_V.size()) +
                                           " values per field and the device's mesh " + std::to_string(m.node_count()) +
                                           " nodes: they do not belong together"));
        }
        return &p;
    };

    ViewOptions view;
    view.contacts = o.contacts;
    view.mesh_edges = o.edges;
    if (o.clip) {
        if (m.dimension() != 3) return std::unexpected(invalid("a clipping box needs a 3D device"));
        if (o.kind != ViewKind::solid && o.kind != ViewKind::field) {
            return std::unexpected(invalid("a clipping box applies to the solid and field views"));
        }
        Box box;
        for (std::size_t k = 0; k < 3; ++k) {
            box.lo[k] = o.clip->lo_um[k] * 1e-4;
            box.hi[k] = o.clip->hi_um[k] * 1e-4;
        }
        view.clip = box;
    }
    // A cut line (Unit 26b, XD5): plotted along by the profile and bands views, drawn on 2D and 3D device views.
    std::optional<CutLine> cut;
    if (o.cut) {
        const bool device_view = o.kind == ViewKind::solid || o.kind == ViewKind::mesh || o.kind == ViewKind::field ||
                                 o.kind == ViewKind::slice || o.kind == ViewKind::iso || o.kind == ViewKind::quality ||
                                 o.kind == ViewKind::difference;
        if (o.kind != ViewKind::profile && o.kind != ViewKind::bands && !(device_view && m.dimension() >= 2)) {
            return std::unexpected(invalid("a cut line applies to the profile and bands views and to 2D and 3D device views"));
        }
        const auto dim = static_cast<std::size_t>(m.dimension());
        if (o.cut->from_um.size() != dim || o.cut->to_um.size() != dim) {
            return std::unexpected(invalid("a cut line needs its end points in the device's dimension"));
        }
        Point from{}, to{};
        for (std::size_t k = 0; k < dim; ++k) {
            from[k] = o.cut->from_um[k] * 1e-4;
            to[k] = o.cut->to_um[k] * 1e-4;
        }
        auto c = cut_line(geometry, from, to);
        if (!c) return std::unexpected(std::move(c.error()));
        cut = std::move(*c);
    }
    view.transparent_insulators = !o.opaque;
    view.projection = o.perspective ? Projection::perspective : Projection::orthographic;
    if (o.azimuth) view.azimuth_deg = *o.azimuth;
    if (o.elevation) view.elevation_deg = *o.elevation;

    const auto field = find_field(o.field);
    if (!field) return std::unexpected(invalid("unknown field " + o.field));
    const auto scale_kind = [&]() -> std::optional<ScaleKind> {
        if (o.scale.empty()) return field_info(*field).default_scale;
        if (o.scale == "linear") return ScaleKind::linear;
        if (o.scale == "log") return ScaleKind::log10;
        if (o.scale == "signed-log") return ScaleKind::signed_log10;
        return std::nullopt;
    }();
    if (!scale_kind) return std::unexpected(invalid("unknown scale " + o.scale));
    const auto colormap = find_colormap(o.colormap);
    if (!colormap) return std::unexpected(invalid("unknown colormap " + o.colormap));
    // Arrows and streamlines (Unit 26c, YD1), with `nitcad plot`'s messages.
    if (!o.streamlines.empty() && !vector_field(o.streamlines)) return std::unexpected(invalid("--streamlines is E, J, Jn or Jp"));
    if (!o.vectors.empty() && !vector_field(o.vectors)) return std::unexpected(invalid("--vectors is E, J, Jn or Jp"));
    if (o.arrows > 100000) return std::unexpected(invalid("--arrows needs a count from 1 to 100000"));
    if (o.lines > 10000) return std::unexpected(invalid("--lines needs a count from 1 to 10000"));
    const std::size_t lines = o.lines > 0 ? o.lines : 24;
    // 3D arrows or streamlines over their own magnitude: on the solid device, or the lines in its outline,
    // colored by it.
    const bool colored_arrows = o.kind == ViewKind::field && m.dimension() == 3 && !o.vectors.empty() &&
                                o.streamlines.empty() && o.field == o.vectors;
    const bool colored_lines = o.kind == ViewKind::field && m.dimension() == 3 && !o.streamlines.empty() &&
                               o.field == o.streamlines;

    const auto values_for = [&]() -> std::expected<std::vector<double>, base::Error> {
        if (*field == Field::net_doping || *field == Field::donors || *field == Field::acceptors) {
            return field_values(device, results::BiasPoint{}, *field);
        }
        auto p = pick();
        if (!p) return std::unexpected(std::move(p.error()));
        auto v = field_values(device, geometry, **p, *field);
        if (!v) return std::unexpected(std::move(v.error()));
        // A vector magnitude on a log scale: without a range, the ten decades below its largest value.
        if (is_vector_magnitude(*field) && *scale_kind == ScaleKind::log10 && !o.range) {
            double top = 0.0;
            for (const double x : *v) if (std::isfinite(x)) top = std::max(top, x);
            if (top > 0.0) {
                const double floor = std::pow(10.0, std::ceil(std::log10(top)) - 10.0);
                for (double& x : *v) if (std::isfinite(x) && x < floor) x = floor;
            }
        }
        return *v;
    };
    const auto field_options = [&](const std::vector<double>& values) -> std::expected<FieldOptions, base::Error> {
        FieldOptions f;
        f.view = view;
        f.map = *colormap;
        f.contours = o.contours;
        f.title = std::string(field_info(*field).title);
        const bool all_zero = is_vector_magnitude(*field) &&
                              std::ranges::none_of(values, [](double x) { return std::isfinite(x) && x != 0.0; });
        auto s = o.range    ? ColorScale::create(*scale_kind, o.range->first, o.range->second, 1e10)
                 : all_zero ? ColorScale::create(ScaleKind::linear, 0.0, 1.0, 1e10)
                            : ColorScale::fit(*scale_kind, values, 1e10);
        if (!s) return std::unexpected(std::move(s.error()));
        f.scale = *s;
        return f;
    };
    std::vector<Plane> planes;
    for (const auto& s : o.slices) {
        Plane p;
        p.normal[static_cast<std::size_t>(s.axis)] = 1.0;
        p.offset = s.at_um * 1e-4;
        planes.push_back(p);
    }
    // Streamlines, then arrows, over a device view (Unit 26c, YD1), as `nitcad plot` adds them.
    const auto overlays = [&](FieldScene* scene) -> std::expected<void, base::Error> {
        if (!o.streamlines.empty()) {
            if (!scene || m.dimension() < 2) {
                return std::unexpected(invalid("streamlines need a 2D or 3D device and a field, slice or iso view"));
            }
            const VectorField vf = *vector_field(o.streamlines);
            auto p = pick();
            if (!p) return std::unexpected(std::move(p.error()));
            auto ev = element_vectors(device, geometry, **p, vf);
            if (!ev) return std::unexpected(std::move(ev.error()));
            const std::vector<Point> nodal = streamline_field(device, geometry, *ev, vf);
            const Locator locator(geometry);
            const auto contacts = device.contacts();
            std::vector<Point> seeds;
            if (!o.seeds.empty()) {
                const std::string_view s = o.seeds;
                if (s.size() > 2 && s[1] == '=' && std::string_view("xyz").find(s[0]) != std::string_view::npos) {
                    const auto at = number(s.substr(2));
                    if (!at) return std::unexpected(invalid("--seeds needs a segment or x=<um>, y=<um> or z=<um>"));
                    if (m.dimension() != 3) return std::unexpected(invalid("a seed plane needs a 3D device"));
                    Plane plane;
                    plane.normal[static_cast<std::size_t>(s[0] - 'x')] = 1.0;
                    plane.offset = *at * 1e-4;
                    seeds = plane_seeds(geometry, plane, lines);
                    if (seeds.empty()) return std::unexpected(invalid("the seed plane misses the device"));
                } else {
                    const auto ab = segment_um(s, m.dimension());
                    if (!ab) {
                        return std::unexpected(
                            invalid("--seeds needs <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] in the device's dimension"));
                    }
                    seeds = segment_seeds(ab->first, ab->second, lines);
                }
            } else {
                // From a contact: the one named, else for currents the one carrying the most, for E the first
                // gate or electrode, else the most biased.
                std::size_t c = 0;
                if (!o.from.empty()) {
                    const device::Contact* named = device.find_contact(o.from);
                    if (!named) return std::unexpected(invalid("there is no contact " + o.from));
                    c = static_cast<std::size_t>(named - contacts.data());
                } else if (vf != VectorField::electric) {
                    double best = -1.0;
                    for (std::size_t k = 0; k < contacts.size(); ++k) {
                        double sum = 0.0;
                        for (const double i : contact_node_currents(device, **p, k)) sum += i;
                        if (std::abs(sum) > best) {
                            best = std::abs(sum);
                            c = k;
                        }
                    }
                } else {
                    const auto gate = std::ranges::find_if(contacts, [](const device::Contact& k) {
                        return k.kind != device::ContactKind::ohmic;
                    });
                    if (gate != contacts.end()) {
                        c = static_cast<std::size_t>(gate - contacts.begin());
                    } else {
                        double best = -1.0;
                        for (std::size_t k = 0; k < contacts.size() && k < (*p)->bias_V.size(); ++k) {
                            if (std::abs((*p)->bias_V[k]) > best) {
                                best = std::abs((*p)->bias_V[k]);
                                c = k;
                            }
                        }
                    }
                }
                const std::vector<double> currents =
                    vf == VectorField::electric ? std::vector<double>{} : contact_node_currents(device, **p, c);
                seeds = contact_seeds(device, locator, c, lines, currents);
            }
            std::vector<Streamline> traced;
            for (const Point& seed : seeds) traced.push_back(trace_streamline(locator, nodal, seed));
            if (colored_lines && scene->colorbar) {
                add_streamlines(*scene, traced, *scene->colorbar);
            } else {
                add_streamlines(*scene, traced, Rgba{0.05f, 0.05f, 0.05f, 0.85f});
            }
        }
        if (!o.vectors.empty()) {
            if (!scene || m.dimension() < 2) {
                return std::unexpected(invalid("arrows need a 2D or 3D device and a field, slice, iso or solid view"));
            }
            const VectorField vf = *vector_field(o.vectors);
            auto p = pick();
            if (!p) return std::unexpected(std::move(p.error()));
            auto ev = element_vectors(device, geometry, **p, vf);
            if (!ev) return std::unexpected(std::move(ev.error()));
            const std::vector<Point> nodal = streamline_field(device, geometry, *ev, vf);
            const Locator locator(geometry);
            const std::size_t count = o.arrows > 0 ? o.arrows : (m.dimension() == 2 ? 300 : 400);
            // Currents: no arrows in the insulators, which carry none.
            std::vector<std::uint8_t> insulator;
            if (vf != VectorField::electric) {
                for (std::size_t i = 0; i < m.node_count(); ++i) {
                    insulator.push_back(device.is_insulator(static_cast<mesh::NodeId>(i)) ? 1 : 0);
                }
            }
            double spacing = 0.0;
            std::vector<Arrow> arrows;
            if (m.dimension() == 2) {
                arrows = grid_arrows(locator, nodal, count, spacing, insulator);
            } else {
                std::vector<Polygon> surface;
                if (o.kind == ViewKind::slice) {
                    for (const Plane& plane : planes) {
                        for (Polygon& q : geometry.slice(plane)) surface.push_back(std::move(q));
                    }
                } else {
                    surface = exposed_boundary(device, geometry);  // the contacts' slabs cover theirs
                }
                arrows = surface_arrows(locator, nodal, surface, count, spacing, insulator);
            }
            if (colored_arrows && scene->colorbar) {
                add_arrows(*scene, arrows, 0.8 * spacing, *scene->colorbar);
            } else {
                add_arrows(*scene, arrows, 0.8 * spacing, Rgba{0.05f, 0.05f, 0.05f, 1.0f});
            }
        }
        return {};
    };
    const auto wrap = [&](auto&& s) -> std::expected<Scene, base::Error> {
        if (!s) return std::unexpected(std::move(s.error()));
        Scene scene(std::move(*s));
        auto* f = std::get_if<FieldScene>(&scene);
        if (auto drawn = overlays(f); !drawn) return std::unexpected(std::move(drawn.error()));
        if (f) {
            if (cut) add_cut_line(*f, *cut);
            apply_camera(*f, o);
        }
        return scene;
    };

    switch (o.kind) {
    case ViewKind::solid:
    case ViewKind::mesh:
        return wrap(o.kind == ViewKind::solid ? solid_view(device, geometry, view) : mesh_view(device, geometry, view));
    case ViewKind::field:
    case ViewKind::slice: {
        auto v = values_for();
        if (!v) return std::unexpected(std::move(v.error()));
        auto f = field_options(*v);
        if (!f) return std::unexpected(std::move(f.error()));
        if (colored_arrows || colored_lines) {
            auto s = colored_arrows ? solid_view(device, geometry, view) : outline_view(device, geometry, view);
            if (s) {
                s->title = f->title;
                s->colorbar = Colorbar{f->map, *f->scale, f->title};
            }
            return wrap(std::move(s));
        }
        if (m.dimension() == 1) {
            return wrap(profile_plot(m, *v, std::string(field_info(*field).title),
                                     *scale_kind == ScaleKind::log10 ? AxisScale::log10 : AxisScale::linear));
        }
        if (o.kind == ViewKind::slice) {
            if (o.slices.empty()) return std::unexpected(invalid("--view slice needs --slice"));
            return wrap(slice_view(device, geometry, *v, planes, *f));
        }
        return wrap(field_view(device, geometry, *v, *f));
    }
    case ViewKind::iso: {
        if (o.levels.empty()) return std::unexpected(invalid("--view iso needs --levels"));
        if (m.dimension() == 1) return std::unexpected(invalid("isosurfaces need a 2D or 3D device (a 1D field is a profile)"));
        auto v = values_for();
        if (!v) return std::unexpected(std::move(v.error()));
        auto f = field_options(*v);
        if (!f) return std::unexpected(std::move(f.error()));
        return wrap(iso_view(device, geometry, *v, o.levels, *f));
    }
    case ViewKind::profile:
    case ViewKind::bands: {
        if (m.dimension() != 1 && !cut) {
            return std::unexpected(invalid("--view " + std::string(to_string(o.kind)) + " in 2D and 3D needs --cut"));
        }
        if (o.kind == ViewKind::bands) {
            auto p = pick();
            if (!p) return std::unexpected(std::move(p.error()));
            return wrap(cut ? cut_band_plot(m, *cut, (*p)->bands) : band_plot(m, (*p)->bands));
        }
        // A profile of the doping needs no results.
        auto v = values_for();
        if (!v) return std::unexpected(std::move(v.error()));
        const AxisScale y = o.log ? AxisScale::log10 : AxisScale::linear;
        std::string title(field_info(*field).title);
        return wrap(cut ? cut_plot(m, *cut, *v, std::move(title), y) : profile_plot(m, *v, std::move(title), y));
    }
    case ViewKind::iv:
    case ViewKind::cv:
    case ViewKind::current: {
        const auto contacts = device.contacts();
        std::vector<SeriesInput> series;
        std::string x_label, y_label;
        // The contacts are chosen on run A and found by name in run B (Unit 26b, XD6).
        const auto in_b = [&](const std::string& name) -> std::expected<std::size_t, base::Error> {
            const std::size_t k = contact_index(*b->device, name, b->device->contacts().size());
            if (k == b->device->contacts().size()) return std::unexpected(invalid("run B has no contact " + name));
            return k;
        };
        const auto labelled = [&](std::string name, const char* run) { return b ? name + " (" + run + ")" : name; };
        if (o.kind == ViewKind::iv) {
            if (!w.sweep) return std::unexpected(invalid("an I-V plot needs a sweep's results"));
            std::vector<std::vector<double>> biases;
            for (const auto& p : w.sweep->points) biases.push_back(p.bias_V);
            const std::size_t swept = contact_index(device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(device, o.sense, swept);
            auto c = analysis::current_curve(*w.sweep, swept, sense);
            if (!c) return std::unexpected(std::move(c.error()));
            series.push_back({labelled("I(" + contacts[sense].name + ")", "A"), std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
            if (b) {
                if (!b->sweep) return std::unexpected(invalid("run B holds no sweep to compare"));
                auto bs = in_b(contacts[swept].name);
                if (!bs) return std::unexpected(std::move(bs.error()));
                auto bn = in_b(contacts[sense].name);
                if (!bn) return std::unexpected(std::move(bn.error()));
                auto cb = analysis::current_curve(*b->sweep, *bs, *bn);
                if (!cb) return std::unexpected(std::move(cb.error()));
                series.push_back({labelled("I(" + contacts[sense].name + ")", "B"), std::move(*cb)});
            }
        } else if (o.kind == ViewKind::cv) {
            if (!w.small_signal) return std::unexpected(invalid("a C-V plot needs a small-signal run's results"));
            std::vector<std::vector<double>> biases;
            for (const auto& p : w.small_signal->points) biases.push_back(p.dc.bias_V);
            const std::size_t swept = contact_index(device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(device, o.sense, swept);
            const std::size_t f = o.frequency.value_or(0);
            if (f >= w.small_signal->frequency_Hz.size()) return std::unexpected(invalid("the frequency is out of range"));
            auto c = analysis::capacitance_curve(*w.small_signal, f, sense, swept, swept);
            if (!c) return std::unexpected(std::move(c.error()));
            const std::string name = "C at " + std::to_string(w.small_signal->frequency_Hz[f]) + " Hz";
            series.push_back({labelled(name, "A"), std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = m.dimension() == 1 ? "C (F/cm\xC2\xB2)" : m.dimension() == 2 ? "C (F/cm)" : "C (F)";
            if (b) {
                if (!b->small_signal || f >= b->small_signal->frequency_Hz.size() ||
                    b->small_signal->frequency_Hz[f] != w.small_signal->frequency_Hz[f]) {
                    return std::unexpected(invalid("run B holds no small-signal run at this frequency to compare"));
                }
                auto bs = in_b(contacts[swept].name);
                if (!bs) return std::unexpected(std::move(bs.error()));
                auto bn = in_b(contacts[sense].name);
                if (!bn) return std::unexpected(std::move(bn.error()));
                auto cb = analysis::capacitance_curve(*b->small_signal, f, *bn, *bs, *bs);
                if (!cb) return std::unexpected(std::move(cb.error()));
                series.push_back({labelled(name, "B"), std::move(*cb)});
            }
        } else {
            if (!w.transient) return std::unexpected(invalid("a current plot needs a transient run's results"));
            if (b && !b->transient) return std::unexpected(invalid("run B holds no transient run to compare"));
            for (std::size_t k = 0; k < contacts.size(); ++k) {
                if (!o.contact.empty() && contacts[k].name != o.contact) continue;
                const auto curve_of = [](const results::Transient& t, std::size_t c) {
                    std::vector<double> time, i;
                    for (const auto& p : t.points) {
                        time.push_back(p.time_s);
                        i.push_back(p.terminal_current[c]);
                    }
                    return analysis::make_curve(std::move(time), std::move(i));
                };
                auto c = curve_of(*w.transient, k);
                if (!c) return std::unexpected(std::move(c.error()));
                series.push_back({labelled(contacts[k].name, "A"), std::move(*c)});
                if (b) {
                    auto bk = in_b(contacts[k].name);
                    if (!bk) return std::unexpected(std::move(bk.error()));
                    auto cb = curve_of(*b->transient, *bk);
                    if (!cb) return std::unexpected(std::move(cb.error()));
                    series.push_back({labelled(contacts[k].name, "B"), std::move(*cb)});
                }
            }
            x_label = "t (s)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
        }
        return wrap(curve_plot(series, x_label, y_label, AxisScale::linear, o.log ? AxisScale::log10 : AxisScale::linear));
    }
    case ViewKind::quality:
    case ViewKind::quality_histogram: {
        auto q = element_quality(geometry);
        if (!q) return std::unexpected(std::move(q.error()));
        if (o.kind == ViewKind::quality_histogram) return wrap(quality_plot(*q, m.dimension()));
        FieldOptions f;
        f.view = view;
        f.map = *colormap;
        if (o.range) {
            auto s = ColorScale::create(ScaleKind::linear, o.range->first, o.range->second);
            if (!s) return std::unexpected(std::move(s.error()));
            f.scale = *s;
        }
        return wrap(quality_view(device, geometry, *q, planes, f));
    }
    case ViewKind::difference: {
        if (!b || !b->device) return std::unexpected(invalid("a difference map needs run B (--compare)"));
        if (!comparable(w, *b)) return std::unexpected(invalid("run B was computed on another mesh: a difference needs one mesh"));
        const bool doping = *field == Field::net_doping || *field == Field::donors || *field == Field::acceptors;
        const results::BiasPoint none;
        const results::BiasPoint* pa = &none;
        const results::BiasPoint* pb = &none;
        if (!doping) {
            auto p = pick();
            if (!p) return std::unexpected(std::move(p.error()));
            pa = *p;
            const auto k = matching_point(*pa, b->points);
            if (!k) return std::unexpected(invalid("run B has no bias point at run A's biases"));
            pb = &b->points[*k];
        }
        auto d = difference_field(device, *b->device, geometry, *pa, *pb, *field);
        if (!d) return std::unexpected(std::move(d.error()));
        FieldOptions f;
        f.view = view;
        f.map = o.colormap == "turbo" ? Colormap::coolwarm : *colormap;  // XD6: cool-warm unless one is chosen
        f.contours = o.contours;
        f.title = d->title;
        if (o.range) {
            auto s = ColorScale::create(ScaleKind::linear, o.range->first, o.range->second);
            if (!s) return std::unexpected(std::move(s.error()));
            f.scale = *s;
        } else {
            f.scale = difference_scale(d->values);
        }
        if (m.dimension() == 1) return wrap(profile_plot(m, d->values, d->title));
        if (m.dimension() == 3 && !planes.empty()) return wrap(slice_view(device, geometry, d->values, planes, f));
        return wrap(field_view(device, geometry, d->values, f));
    }
    case ViewKind::exploded: {
        // Unit 26c (YD2): `nitcad plot --view exploded [--axis] [--gap] [--field]`.
        if (m.dimension() != 3) return std::unexpected(invalid("the exploded view needs a 3D device"));
        if (o.axis < 0 || o.axis > 2) return std::unexpected(invalid("--axis is x, y or z"));
        if (o.gap_um && !(*o.gap_um >= 0.0)) return std::unexpected(invalid("--gap needs a length in um, at least 0"));
        const double gap = o.gap_um ? *o.gap_um * 1e-4 : default_gap(device, o.axis);
        if (!o.exploded_field) {
            FieldOptions f;
            f.view = view;
            return wrap(exploded_view(device, geometry, o.axis, gap, {}, f));
        }
        auto v = values_for();
        if (!v) return std::unexpected(std::move(v.error()));
        auto f = field_options(*v);
        if (!f) return std::unexpected(std::move(f.error()));
        return wrap(exploded_view(device, geometry, o.axis, gap, *v, *f));
    }
    }
    return std::unexpected(invalid("unknown view"));
}

std::vector<ViewKind> available_views(const Workspace& w, const Workspace* b) {
    std::vector<ViewKind> out;
    if (!w.device) return out;
    const int dim = w.device->mesh().dimension();
    const bool elements = w.geometry && w.geometry->element_count() > 0;
    if (dim >= 2 && elements) {
        out.push_back(ViewKind::solid);
        out.push_back(ViewKind::mesh);
    }
    out.push_back(dim == 1 ? ViewKind::profile : ViewKind::field);
    if (dim == 3 && elements) out.push_back(ViewKind::slice);
    if (dim == 3 && elements) out.push_back(ViewKind::exploded);  // Unit 26c (YD2)
    if (dim >= 2 && elements && w.has_results()) out.push_back(ViewKind::iso);
    if (dim == 1 && w.has_results()) out.push_back(ViewKind::bands);
    // Along a cut line (Unit 26b, XD5).
    if (dim >= 2 && elements) out.push_back(ViewKind::profile);
    if (dim >= 2 && elements && w.has_results()) out.push_back(ViewKind::bands);
    if (w.sweep && w.sweep->points.size() >= 2) out.push_back(ViewKind::iv);
    if (w.small_signal && w.small_signal->points.size() >= 2) out.push_back(ViewKind::cv);
    if (w.transient && w.transient->points.size() >= 2) out.push_back(ViewKind::current);
    // Unit 26b: mesh quality (XD7) and the difference with run B (XD6).
    if (dim >= 2 && elements) {
        out.push_back(ViewKind::quality);
        out.push_back(ViewKind::quality_histogram);
    }
    if (b && comparable(w, *b)) out.push_back(ViewKind::difference);
    return out;
}

std::vector<std::string> available_fields(const Workspace& w) {
    std::vector<std::string> out;
    if (!w.device) return out;
    const bool results = w.has_results();
    const results::BiasPoint* last = results ? &w.points.back() : nullptr;
    const bool elements = w.geometry && w.geometry->element_count() > 0;
    for (int k = 0; k <= static_cast<int>(Field::last); ++k) {
        const Field f = static_cast<Field>(k);
        const bool doping = f == Field::net_doping || f == Field::donors || f == Field::acceptors;
        if (!doping && !results) continue;
        if (is_vector_magnitude(f) && !elements) continue;
        if (!doping && last) {
            // The rates only when the result stores them: ask for the values and keep what works.
            if (f >= Field::recombination && f <= Field::btbt_generation &&
                !field_values(*w.device, *last, f).has_value()) {
                continue;
            }
        }
        out.emplace_back(field_info(f).key);
    }
    return out;
}

namespace {

// render::ViewOptions' defaults, the angles a setting starts from.
constexpr float default_azimuth = 35.0f, default_elevation = 25.0f;

}  // namespace

void orbit(ViewSettings& s, float azimuth_deg, float elevation_deg) {
    float az = s.azimuth.value_or(default_azimuth) + azimuth_deg;
    az = std::fmod(az, 360.0f);
    if (az > 180.0f) az -= 360.0f;
    if (az <= -180.0f) az += 360.0f;
    s.azimuth = az;
    s.elevation = std::clamp(s.elevation.value_or(default_elevation) + elevation_deg, -89.0f, 89.0f);
}

void zoom_by(ViewSettings& s, float factor) {
    if (!(factor > 0.0f) || !std::isfinite(factor)) return;
    s.zoom = std::clamp(s.zoom * factor, 1.0f / 16.0f, 1024.0f);
}

void pan_by(ViewSettings& s, float dx, float dy) {
    if (!std::isfinite(dx) || !std::isfinite(dy)) return;
    s.pan_x += dx;
    s.pan_y += dy;
}

void fit(ViewSettings& s) {
    s.zoom = 1.0f;
    s.pan_x = s.pan_y = 0.0f;
}

void apply_camera(FieldScene& scene, const ViewSettings& s) {
    if (scene.dimension == 3) {
        scene.camera.azimuth_deg = s.azimuth.value_or(default_azimuth);
        scene.camera.elevation_deg = s.elevation.value_or(default_elevation);
        scene.camera.projection = s.perspective ? Projection::perspective : Projection::orthographic;
    }
    scene.camera.zoom = s.zoom;
    scene.camera.pan = {s.pan_x, s.pan_y};
}

std::vector<std::string> plot_arguments(const ViewSettings& s, const Workspace* b) {
    std::vector<std::string> a{"--view", std::string(to_string(s.kind))};
    const bool fieldish = s.kind == ViewKind::field || s.kind == ViewKind::slice || s.kind == ViewKind::iso ||
                          s.kind == ViewKind::profile || s.kind == ViewKind::difference ||
                          (s.kind == ViewKind::exploded && s.exploded_field);
    const bool uses_b = s.kind == ViewKind::difference || s.kind == ViewKind::iv || s.kind == ViewKind::cv ||
                        s.kind == ViewKind::current;
    if (b && uses_b) a.insert(a.end(), {"--compare", b->source.string()});
    if (fieldish) a.insert(a.end(), {"--field", s.field});
    if (s.point) a.insert(a.end(), {"--point", std::to_string(*s.point)});
    if (!s.scale.empty() && fieldish) a.insert(a.end(), {"--scale", s.scale});
    if (s.colormap != "turbo") a.insert(a.end(), {"--colormap", s.colormap});
    if (s.range) a.insert(a.end(), {"--range", std::format("{},{}", s.range->first, s.range->second)});
    if (s.contours) a.emplace_back("--contours");
    if (s.edges) a.emplace_back("--edges");
    if (s.opaque) a.emplace_back("--opaque");
    if (s.perspective) a.emplace_back("--perspective");
    if (s.azimuth) a.insert(a.end(), {"--azimuth", std::format("{}", *s.azimuth)});
    if (s.elevation) a.insert(a.end(), {"--elevation", std::format("{}", *s.elevation)});
    for (const auto& p : s.slices) {
        a.insert(a.end(), {"--slice", std::string(1, static_cast<char>('x' + p.axis)) + "=" + std::format("{}", p.at_um)});
    }
    if (s.clip) {
        const auto& c = *s.clip;
        a.insert(a.end(), {"--clip", std::format("{},{},{}:{},{},{}", c.lo_um[0], c.lo_um[1], c.lo_um[2], c.hi_um[0],
                                                 c.hi_um[1], c.hi_um[2])});
    }
    if (!s.levels.empty()) {
        std::string l;
        for (const double v : s.levels) l += (l.empty() ? "" : ",") + std::format("{}", v);
        a.insert(a.end(), {"--levels", l});
    }
    if (!s.contact.empty()) a.insert(a.end(), {"--contact", s.contact});
    if (!s.sense.empty()) a.insert(a.end(), {"--sense", s.sense});
    if (s.frequency) a.insert(a.end(), {"--frequency", std::to_string(*s.frequency)});
    if (s.log) a.emplace_back("--log");
    if (s.cut) {
        const auto end = [](const std::vector<double>& v) {
            std::string t;
            for (const double x : v) t += (t.empty() ? "" : ",") + std::format("{}", x);
            return t;
        };
        a.insert(a.end(), {"--cut", end(s.cut->from_um) + ":" + end(s.cut->to_um)});
    }
    // Unit 26c (YD1, YD2): arrows and streamlines on the field views (takes_vectors).
    if (takes_vectors(s.kind) && !s.streamlines.empty()) {
        a.insert(a.end(), {"--streamlines", s.streamlines});
        if (s.lines > 0) a.insert(a.end(), {"--lines", std::to_string(s.lines)});
        if (!s.from.empty()) a.insert(a.end(), {"--from", s.from});
        if (!s.seeds.empty()) a.insert(a.end(), {"--seeds", s.seeds});
    }
    if (takes_vectors(s.kind) && !s.vectors.empty()) {
        a.insert(a.end(), {"--vectors", s.vectors});
        if (s.arrows > 0) a.insert(a.end(), {"--arrows", std::to_string(s.arrows)});
    }
    if (s.kind == ViewKind::exploded) {
        if (s.axis != 1) a.insert(a.end(), {"--axis", std::string(1, static_cast<char>('x' + s.axis))});
        if (s.gap_um) a.insert(a.end(), {"--gap", std::format("{}", *s.gap_um)});
    }
    return a;
}

}  // namespace NiTCAD::app
