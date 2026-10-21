// nitcad plot (Unit 25a, DECISIONS.md VR6, VG5): views and curves to PNG images.
#include "plot.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::cli {

namespace {

namespace fs = std::filesystem;
using namespace NiTCAD::render;

constexpr const char* plot_usage =
    "usage: nitcad plot <deck | result file> [options] --out <image.png>\n"
    "  input   a deck (its device, mesh and the result files its run wrote), or a .nitsweep,\n"
    "          .nittran or .nitac file with --device <file.nitdev> [--mesh <file.nitmesh>]\n"
    "  view    --view solid|mesh|field|slice|iv|cv|current|bands|profile\n"
    "          (default: field with --field, else solid)\n"
    "  field   --field potential|n|p|Ec|Ev|Efn|Efp|T|doping|ND|NA  --point <k> (default: last)\n"
    "          --scale linear|log|signed-log  --range <lo>,<hi>\n"
    "          --colormap turbo|viridis|cividis|coolwarm|soft  --contours  --edges  --opaque\n"
    "  3D      --slice x=<um>|y=<um>|z=<um> (repeatable)  --azimuth <deg>  --elevation <deg>\n"
    "          --perspective\n"
    "  curves  --contact <name> (swept)  --sense <name>  --frequency <k>  --log\n"
    "  image   --size <w>x<h> (1280x960)  --adapter auto|gpu|warp\n";

int usage(std::string_view why) {
    std::fprintf(stderr, "nitcad plot: %.*s\n%s", static_cast<int>(why.size()), why.data(), plot_usage);
    return 64;
}

int fail(std::string_view message) {
    std::fprintf(stderr, "nitcad plot: %.*s\n", static_cast<int>(message.size()), message.data());
    return 1;
}

std::optional<double> number(std::string_view s) {
    double v = 0.0;
    const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{} || p != s.data() + s.size() || !std::isfinite(v)) return std::nullopt;
    return v;
}

struct Options {
    fs::path input, device, mesh, out;
    std::string view, field, scale, colormap = "turbo", contact, sense, adapter = "auto";
    std::optional<std::size_t> point, frequency;
    std::optional<std::pair<double, double>> range;
    std::vector<Plane> slices;
    std::optional<float> azimuth, elevation;
    bool contours = false, edges = false, opaque = false, perspective = false, log = false;
    int width = 1280, height = 960;
};

// The results of a run, as the views read them.
struct Loaded {
    std::optional<device::Device> device;  // when not from a deck
    std::optional<mesh::GeneratedMesh> generated;
    std::optional<results::Sweep> sweep;
    std::optional<results::Transient> transient;
    std::optional<results::SmallSignal> small_signal;
    std::optional<results::BiasPoint> equilibrium;
};

// The bias points of whatever was loaded: a sweep's, a transient's snapshots, a small-signal run's
// operating points, or the equilibrium.
std::vector<results::BiasPoint> points_of(const Loaded& l) {
    std::vector<results::BiasPoint> out;
    if (l.sweep) out = l.sweep->points;
    if (l.transient) {
        for (const auto& s : l.transient->snapshots) {
            results::BiasPoint p;
            p.fields = s.fields;
            p.bands = s.bands;
            out.push_back(std::move(p));
        }
    }
    if (l.small_signal) {
        for (const auto& p : l.small_signal->points) out.push_back(p.dc);
    }
    if (l.equilibrium) out.push_back(*l.equilibrium);
    return out;
}

std::expected<void, std::string> load_result(const fs::path& path, Loaded& l) {
    const std::string ext = path.extension().string();
    if (!fs::exists(path)) return std::unexpected("no result file " + path.string() + " (run the deck first)");
    if (ext == ".nitsweep") {
        auto s = io::load_sweep(path);
        if (!s) return std::unexpected(s.error().message);
        l.sweep = std::move(*s);
    } else if (ext == ".nittran") {
        auto t = io::load_transient(path);
        if (!t) return std::unexpected(t.error().message);
        l.transient = std::move(*t);
    } else if (ext == ".nitac") {
        auto a = io::load_small_signal(path);
        if (!a) return std::unexpected(a.error().message);
        l.small_signal = std::move(*a);
    } else {
        return std::unexpected("a result file is .nitsweep, .nittran or .nitac: " + path.string());
    }
    return {};
}

std::string current_unit(int dimension) {
    return dimension == 1 ? "A/cm\xC2\xB2" : dimension == 2 ? "A/cm" : "A";
}

std::size_t contact_index(const device::Device& d, const std::string& name, std::size_t fallback) {
    for (std::size_t k = 0; k < d.contacts().size(); ++k) {
        if (d.contacts()[k].name == name) return k;
    }
    return fallback;
}

// The contact whose bias varies most over the points (the swept one).
std::size_t swept_contact(const std::vector<std::vector<double>>& biases) {
    std::size_t best = 0;
    double spread = -1.0;
    for (std::size_t k = 0; !biases.empty() && k < biases.front().size(); ++k) {
        double lo = 1e300, hi = -1e300;
        for (const auto& b : biases) {
            lo = std::min(lo, b[k]);
            hi = std::max(hi, b[k]);
        }
        if (hi - lo > spread) {
            spread = hi - lo;
            best = k;
        }
    }
    return best;
}

}  // namespace

int plot(int argc, char** argv) {
    Options o;
    if (argc < 3) return usage("missing the input");
    o.input = argv[2];
    for (int i = 3; i < argc; ++i) {
        const std::string_view a = argv[i];
        const auto value = [&]() -> std::optional<std::string_view> {
            if (i + 1 >= argc) return std::nullopt;
            return std::string_view(argv[++i]);
        };
        const auto flag_value = [&](auto& target) -> bool {
            const auto v = value();
            if (!v) return false;
            target = std::string(*v);
            return true;
        };
        if (a == "--out") { if (!flag_value(o.out)) return usage("--out needs a path"); }
        else if (a == "--device") { if (!flag_value(o.device)) return usage("--device needs a path"); }
        else if (a == "--mesh") { if (!flag_value(o.mesh)) return usage("--mesh needs a path"); }
        else if (a == "--view") { if (!flag_value(o.view)) return usage("--view needs a name"); }
        else if (a == "--field") { if (!flag_value(o.field)) return usage("--field needs a name"); }
        else if (a == "--scale") { if (!flag_value(o.scale)) return usage("--scale needs a name"); }
        else if (a == "--colormap") { if (!flag_value(o.colormap)) return usage("--colormap needs a name"); }
        else if (a == "--contact") { if (!flag_value(o.contact)) return usage("--contact needs a name"); }
        else if (a == "--sense") { if (!flag_value(o.sense)) return usage("--sense needs a name"); }
        else if (a == "--adapter") { if (!flag_value(o.adapter)) return usage("--adapter needs a name"); }
        else if (a == "--contours") o.contours = true;
        else if (a == "--edges") o.edges = true;
        else if (a == "--opaque") o.opaque = true;
        else if (a == "--perspective") o.perspective = true;
        else if (a == "--log") o.log = true;
        else if (a == "--point" || a == "--frequency") {
            const auto v = value();
            const auto n = v ? number(*v) : std::nullopt;
            if (!n || *n < 0 || *n != std::floor(*n)) return usage(std::string(a) + " needs an index");
            (a == "--point" ? o.point : o.frequency) = static_cast<std::size_t>(*n);
        } else if (a == "--azimuth" || a == "--elevation") {
            const auto v = value();
            const auto n = v ? number(*v) : std::nullopt;
            if (!n) return usage(std::string(a) + " needs degrees");
            (a == "--azimuth" ? o.azimuth : o.elevation) = static_cast<float>(*n);
        } else if (a == "--range") {
            const auto v = value();
            const auto comma = v ? v->find(',') : std::string_view::npos;
            if (comma == std::string_view::npos) return usage("--range needs <lo>,<hi>");
            const auto lo = number(v->substr(0, comma)), hi = number(v->substr(comma + 1));
            if (!lo || !hi) return usage("--range needs <lo>,<hi>");
            o.range = std::pair{*lo, *hi};
        } else if (a == "--slice") {
            const auto v = value();
            if (!v || v->size() < 3 || (*v)[1] != '=' || std::string_view("xyz").find((*v)[0]) == std::string_view::npos) {
                return usage("--slice needs x=<um>, y=<um> or z=<um>");
            }
            const auto at = number(v->substr(2));
            if (!at) return usage("--slice needs x=<um>, y=<um> or z=<um>");
            Plane p;
            p.normal[static_cast<std::size_t>((*v)[0] - 'x')] = 1.0;
            p.offset = *at * 1e-4;
            o.slices.push_back(p);
        } else if (a == "--size") {
            const auto v = value();
            const auto x = v ? v->find('x') : std::string_view::npos;
            const auto w = x != std::string_view::npos ? number(v->substr(0, x)) : std::nullopt;
            const auto h = x != std::string_view::npos ? number(v->substr(x + 1)) : std::nullopt;
            if (!w || !h) return usage("--size needs <w>x<h>");
            o.width = static_cast<int>(*w);
            o.height = static_cast<int>(*h);
        } else {
            return usage("unknown option " + std::string(a));
        }
    }
    if (o.out.empty()) return usage("--out is required");
    if (o.view.empty()) o.view = !o.slices.empty() ? "slice" : o.field.empty() ? "solid" : "field";
    if (o.field.empty()) o.field = "potential";

    // The device, its geometry and the results.
    Loaded l;
    const device::Device* device = nullptr;
    std::optional<deck::Plan> plan;
    const std::string ext = o.input.extension().string();
    if (ext == ".nitsweep" || ext == ".nittran" || ext == ".nitac") {
        if (o.device.empty()) return usage("a result file needs --device <file.nitdev>");
        auto d = io::load_device(o.device);
        if (!d) return fail(d.error().message);
        l.device = std::move(*d);
        device = &*l.device;
        if (!o.mesh.empty()) {
            auto g = io::load_generated_mesh(o.mesh);
            if (!g) return fail(g.error().message);
            if (mesh::graph_digest(g->mesh) != mesh::graph_digest(device->mesh())) {
                return fail("the mesh file " + o.mesh.string() + " is not the device's mesh");
            }
            l.generated = std::move(*g);
        }
        if (auto r = load_result(o.input, l); !r) return fail(r.error());
    } else {
        const auto deck = deck::read(o.input);
        if (!deck) return fail(deck.error().message);
        auto built = deck::build(*deck);
        if (!built) return fail(built.error().message);
        plan = std::move(*built);
        device = &plan->device;
        const bool needs_results = o.view != "solid" && o.view != "mesh" &&
                                   !(o.view == "field" && (o.field == "doping" || o.field == "ND" || o.field == "NA"));
        if (needs_results) {
            const deck::Outputs& out = plan->outputs;
            using deck::Analysis;
            const fs::path file = plan->analysis == Analysis::transient      ? out.transient
                                  : plan->analysis == Analysis::small_signal ? out.small_signal
                                  : plan->analysis == Analysis::equilibrium  ? out.checkpoint
                                                                             : out.sweep;
            if (file.empty()) return fail("the deck writes no result file for this view (add one to [output])");
            if (plan->analysis == Analysis::equilibrium) {
                if (!fs::exists(file)) return fail("no result file " + file.string() + " (run the deck first)");
                auto c = io::load_checkpoint(file, device->mesh());
                if (!c) return fail(c.error().message);
                results::BiasPoint p;
                p.fields = c->fields;
                l.equilibrium = std::move(p);
            } else if (auto r = load_result(file, l); !r) {
                return fail(r.error());
            }
        }
    }
    const mesh::Mesh& m = device->mesh();
    const mesh::GeneratedMesh* generated = plan && plan->generated ? &*plan->generated : l.generated ? &*l.generated : nullptr;
    std::optional<Geometry> geometry;
    if (device->cells()) {
        auto g = Geometry::tensor(m, *device->cells());
        if (!g) return fail(g.error().message);
        geometry = std::move(*g);
    } else if (generated) {
        auto g = Geometry::simplex(m, generated->simplices);
        if (!g) return fail(g.error().message);
        geometry = std::move(*g);
    } else {
        geometry = Geometry::graph(m);
    }

    // The bias point, checked against the mesh.
    const std::vector<results::BiasPoint> points = points_of(l);
    const auto pick = [&]() -> std::expected<const results::BiasPoint*, std::string> {
        if (points.empty()) return std::unexpected(std::string("the results hold no field to plot"));
        const std::size_t k = o.point.value_or(points.size() - 1);
        if (k >= points.size()) {
            return std::unexpected("--point " + std::to_string(k) + " is out of range (" +
                                   std::to_string(points.size()) + " points)");
        }
        const auto& p = points[k];
        if (p.fields.potential_V.size() != m.node_count()) {
            return std::unexpected("the result has " + std::to_string(p.fields.potential_V.size()) +
                                   " values per field and the device's mesh " + std::to_string(m.node_count()) +
                                   " nodes: they do not belong together");
        }
        return &p;
    };

    // The scene.
    ViewOptions view;
    view.contacts = true;
    view.mesh_edges = o.edges;
    view.transparent_insulators = !o.opaque;
    view.projection = o.perspective ? Projection::perspective : Projection::orthographic;
    if (o.azimuth) view.azimuth_deg = *o.azimuth;
    if (o.elevation) view.elevation_deg = *o.elevation;
    std::optional<FieldScene> field_scene;
    std::optional<PlotScene> plot_scene;
    const auto field = find_field(o.field);
    if (!field) return usage("unknown field " + o.field);
    const auto scale_kind = [&]() -> std::optional<ScaleKind> {
        if (o.scale.empty()) return field_info(*field).default_scale;
        if (o.scale == "linear") return ScaleKind::linear;
        if (o.scale == "log") return ScaleKind::log10;
        if (o.scale == "signed-log") return ScaleKind::signed_log10;
        return std::nullopt;
    }();
    if (!scale_kind) return usage("unknown scale " + o.scale);
    const auto colormap = find_colormap(o.colormap);
    if (!colormap) return usage("unknown colormap " + o.colormap);

    const auto values_for = [&]() -> std::expected<std::vector<double>, std::string> {
        if (*field == Field::net_doping || *field == Field::donors || *field == Field::acceptors) {
            auto v = field_values(*device, results::BiasPoint{}, *field);
            if (!v) return std::unexpected(v.error().message);
            return *v;
        }
        auto p = pick();
        if (!p) return std::unexpected(p.error());
        auto v = field_values(*device, **p, *field);
        if (!v) return std::unexpected(v.error().message);
        return *v;
    };
    const auto field_options = [&](const std::vector<double>& values) -> std::expected<FieldOptions, std::string> {
        FieldOptions f;
        f.view = view;
        f.map = *colormap;
        f.contours = o.contours;
        f.title = std::string(field_info(*field).title);
        auto s = o.range ? ColorScale::create(*scale_kind, o.range->first, o.range->second, 1e10)
                         : ColorScale::fit(*scale_kind, values, 1e10);
        if (!s) return std::unexpected(s.error().message);
        f.scale = *s;
        return f;
    };

    if (o.view == "solid" || o.view == "mesh") {
        auto s = o.view == "solid" ? solid_view(*device, *geometry, view) : mesh_view(*device, *geometry, view);
        if (!s) return fail(s.error().message);
        field_scene = std::move(*s);
    } else if (o.view == "field" || o.view == "slice") {
        auto v = values_for();
        if (!v) return fail(v.error());
        auto f = field_options(*v);
        if (!f) return fail(f.error());
        if (m.dimension() == 1) {
            auto s = profile_plot(m, *v, std::string(field_info(*field).title),
                                  *scale_kind == ScaleKind::log10 ? AxisScale::log10 : AxisScale::linear);
            if (!s) return fail(s.error().message);
            plot_scene = std::move(*s);
        } else if (o.view == "slice") {
            if (o.slices.empty()) return usage("--view slice needs --slice");
            auto s = slice_view(*device, *geometry, *v, o.slices, *f);
            if (!s) return fail(s.error().message);
            field_scene = std::move(*s);
        } else {
            auto s = field_view(*device, *geometry, *v, *f);
            if (!s) return fail(s.error().message);
            field_scene = std::move(*s);
        }
    } else if (o.view == "profile" || o.view == "bands") {
        auto p = pick();
        if (!p) return fail(p.error());
        if (o.view == "bands") {
            auto s = band_plot(m, (*p)->bands);
            if (!s) return fail(s.error().message);
            plot_scene = std::move(*s);
        } else {
            auto v = values_for();
            if (!v) return fail(v.error());
            auto s = profile_plot(m, *v, std::string(field_info(*field).title), o.log ? AxisScale::log10 : AxisScale::linear);
            if (!s) return fail(s.error().message);
            plot_scene = std::move(*s);
        }
    } else if (o.view == "iv" || o.view == "cv" || o.view == "current") {
        const auto contacts = device->contacts();
        std::vector<SeriesInput> series;
        std::string x_label, y_label;
        AxisScale x_scale = AxisScale::linear;
        if (o.view == "iv") {
            if (!l.sweep) return fail("an I-V plot needs a sweep's results");
            std::vector<std::vector<double>> biases;
            for (const auto& p : l.sweep->points) biases.push_back(p.bias_V);
            const std::size_t swept = contact_index(*device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(*device, o.sense, swept);
            auto c = analysis::current_curve(*l.sweep, swept, sense);
            if (!c) return fail(c.error().message);
            series.push_back({"I(" + contacts[sense].name + ")", std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
        } else if (o.view == "cv") {
            if (!l.small_signal) return fail("a C-V plot needs a small-signal run's results");
            std::vector<std::vector<double>> biases;
            for (const auto& p : l.small_signal->points) biases.push_back(p.dc.bias_V);
            const std::size_t swept = contact_index(*device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(*device, o.sense, swept);
            const std::size_t f = o.frequency.value_or(0);
            auto c = analysis::capacitance_curve(*l.small_signal, f, sense, swept, swept);
            if (!c) return fail(c.error().message);
            series.push_back({"C at " + std::to_string(l.small_signal->frequency_Hz[f]) + " Hz", std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = m.dimension() == 1 ? "C (F/cm\xC2\xB2)" : m.dimension() == 2 ? "C (F/cm)" : "C (F)";
        } else {
            if (!l.transient) return fail("a current plot needs a transient run's results");
            for (std::size_t k = 0; k < contacts.size(); ++k) {
                if (!o.contact.empty() && contacts[k].name != o.contact) continue;
                std::vector<double> t, i;
                for (const auto& p : l.transient->points) {
                    t.push_back(p.time_s);
                    i.push_back(p.terminal_current[k]);
                }
                auto c = analysis::make_curve(std::move(t), std::move(i));
                if (!c) return fail(c.error().message);
                series.push_back({contacts[k].name, std::move(*c)});
            }
            x_label = "t (s)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
        }
        auto s = curve_plot(series, x_label, y_label, x_scale, o.log ? AxisScale::log10 : AxisScale::linear);
        if (!s) return fail(s.error().message);
        plot_scene = std::move(*s);
    } else {
        return usage("unknown view " + o.view);
    }

    // Draw and write.
    const Adapter adapter = o.adapter == "gpu" ? Adapter::hardware : o.adapter == "warp" ? Adapter::warp : Adapter::automatic;
    if (o.adapter != "gpu" && o.adapter != "warp" && o.adapter != "auto") return usage("unknown adapter " + o.adapter);
    auto renderer = Renderer::create({.adapter = adapter});
    if (!renderer) return fail(renderer.error().message);
    auto image = field_scene ? renderer->draw(*field_scene, o.width, o.height) : renderer->draw(*plot_scene, o.width, o.height);
    if (!image) return fail(image.error().message);
    if (auto saved = save_png(*image, o.out); !saved) return fail(saved.error().message);
    std::printf("wrote %s (%dx%d, %s)\n", o.out.string().c_str(), image->width, image->height,
                renderer->adapter_name().c_str());
    return 0;
}

const char* plot_usage_text() { return plot_usage; }

}  // namespace NiTCAD::cli
