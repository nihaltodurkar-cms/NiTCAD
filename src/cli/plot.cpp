// nitcad plot (Unit 25a, DECISIONS.md VR6, VG5): views and curves to PNG images.
#include "plot.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/render/compare.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/quality.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/render/vectors.hpp"

namespace NiTCAD::cli {

namespace {

namespace fs = std::filesystem;
using namespace NiTCAD::render;

constexpr const char* plot_usage =
    "usage: nitcad plot <deck | result file> [options] --out <image.png>\n"
    "  input   a deck (its device, mesh and the result files its run wrote), or a .nitsweep,\n"
    "          .nittran or .nitac file with --device <file.nitdev> [--mesh <file.nitmesh>]\n"
    "  view    --view solid|mesh|field|slice|iso|exploded|filmstrip|iv|cv|current|bands|profile|\n"
    "          quality|quality-histogram|difference (default: field with --field, else solid)\n"
    "  iso     --view iso --field <f> --levels <v>[,<v>...] (3D surfaces, 2D contours)\n"
    "  field   --field potential|n|p|Ec|Ev|Efn|Efp|T|doping|ND|NA|R|srh|auger|\n"
    "          radiative|surface|Gii|Gbtbt|E|J|Jn|Jp  --point <k> (default: last)\n"
    "          --scale linear|log|signed-log  --range <lo>,<hi>\n"
    "          --colormap turbo|viridis|cividis|coolwarm|soft  --contours  --edges  --opaque\n"
    "  3D      --slice x=<um>|y=<um>|z=<um> (repeatable)  --azimuth <deg>  --elevation <deg>\n"
    "          --perspective\n"
    "  curves  --contact <name> (swept)  --sense <name>  --frequency <k>  --log\n"
    "  arrows  --vectors E|J|Jn|Jp  --arrows <n> (300 in 2D, 400 in 3D)\n"
    "  clip    --clip <x0>,<y0>,<z0>:<x1>,<y1>,<z1> (um): that box removed from a 3D solid or field view\n"
    "  cut     --cut <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] (um): --view profile or bands along that line (needed\n"
    "          in 2D and 3D), or the line drawn on a 2D or 3D solid, mesh, field, slice or iso view\n"
    "  compare --compare <deck | result file> (run B; a result file is read with the device of run A):\n"
    "          --view difference --field <f>: B - A at the bias point both have, on one mesh (mV for\n"
    "          the potential, meV for energies, decades for densities); with --view iv, cv or current\n"
    "          both runs' curves\n"
    "  quality --view quality (each element's minimum angle, 2D, or dihedral angle, 3D; --slice in 3D)\n"
    "          or --view quality-histogram\n"
    "  explode --view exploded [--field <f>] [--axis x|y|z (y)] [--gap <um> (12% of the extent)]\n"
    "  film    --view filmstrip --field <f> (a transient run) [--every <k>] [--columns <c> (4)]\n"
    "          (--size is each frame's, 640x480 by default)\n"
    "          [--frames <dir>] (each frame also as <dir>/frame_000.png, ...)\n"
    "  lines   --streamlines E|J|Jn|Jp  --lines <n> (24)  --from <contact>\n"
    "          --seeds <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] (um) | x=<um>|y=<um>|z=<um> (3D plane)\n"
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

// A segment "<x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>]" in um, as its ends in cm, each end with exactly `dimension`
// coordinates.
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

struct Options {
    fs::path input, device, mesh, out;
    std::string view, field, scale, colormap = "turbo", contact, sense, adapter = "auto";
    std::string streamlines, from, seeds;  // Unit 25b
    std::string cut;  // Unit 26b
    fs::path compare;  // Unit 26b: run B
    std::size_t lines = 24;
    std::string vectors;  // Unit 25c
    std::optional<Box> clip;
    std::optional<double> gap;  // um
    int axis = 1;
    std::size_t every = 0, columns = 0;  // filmstrip
    fs::path frames;
    std::size_t arrows = 0;
    std::optional<std::size_t> point, frequency;
    std::optional<std::pair<double, double>> range;
    std::vector<Plane> slices;
    std::vector<double> levels;
    std::optional<float> azimuth, elevation;
    bool contours = false, edges = false, opaque = false, perspective = false, log = false;
    int width = 1280, height = 960;
    bool size_given = false;
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
            // Unit 25c: with their conduction currents and rates, for |J|, lines, arrows and R.
            p.edge_current_n = s.edge_current_n;
            p.edge_current_p = s.edge_current_p;
            p.rates = s.rates;
            if (s.point < l.transient->points.size()) p.bias_V = l.transient->points[s.point].bias_V;
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

// The results a deck's run wrote, for its kind of analysis.
std::expected<void, std::string> load_outputs(const deck::Plan& plan, Loaded& l) {
    const deck::Outputs& out = plan.outputs;
    using deck::Analysis;
    const fs::path file = plan.analysis == Analysis::transient      ? out.transient
                          : plan.analysis == Analysis::small_signal ? out.small_signal
                          : plan.analysis == Analysis::equilibrium  ? out.checkpoint
                                                                     : out.sweep;
    if (file.empty()) return std::unexpected(std::string("the deck writes no result file for this view (add one to [output])"));
    if (plan.analysis == Analysis::equilibrium) {
        if (!fs::exists(file)) return std::unexpected("no result file " + file.string() + " (run the deck first)");
        auto c = io::load_checkpoint(file, plan.device.mesh());
        if (!c) return std::unexpected(c.error().message);
        results::BiasPoint p;
        p.fields = c->fields;
        l.equilibrium = std::move(p);
        return {};
    }
    return load_result(file, l);
}

// The mesh digest a loaded result names (0 when it names none: an equilibrium or an older file).
std::uint64_t result_digest(const Loaded& l) {
    return l.sweep ? l.sweep->run.mesh_digest : l.transient ? l.transient->run.mesh_digest
           : l.small_signal ? l.small_signal->run.mesh_digest : 0;
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
        else if (a == "--streamlines") { if (!flag_value(o.streamlines)) return usage("--streamlines needs E, J, Jn or Jp"); }
        else if (a == "--vectors") { if (!flag_value(o.vectors)) return usage("--vectors needs E, J, Jn or Jp"); }
        else if (a == "--frames") { if (!flag_value(o.frames)) return usage("--frames needs a directory"); }
        else if (a == "--every" || a == "--columns") {
            const auto v = value();
            const auto n = v ? number(*v) : std::nullopt;
            if (!n || *n < 1 || *n > 1000 || *n != std::floor(*n)) return usage(std::string(a) + " needs a count from 1 to 1000");
            (a == "--every" ? o.every : o.columns) = static_cast<std::size_t>(*n);
        }
        else if (a == "--gap") {
            const auto v = value();
            const auto n = v ? number(*v) : std::nullopt;
            if (!n || *n < 0.0) return usage("--gap needs a length in um, at least 0");
            o.gap = *n;
        } else if (a == "--axis") {
            const auto v = value();
            if (!v || v->size() != 1 || std::string_view("xyz").find((*v)[0]) == std::string_view::npos) {
                return usage("--axis is x, y or z");
            }
            o.axis = (*v)[0] - 'x';
        }
        else if (a == "--clip") {
            // <x0>,<y0>,<z0>:<x1>,<y1>,<z1> in um, lo below hi on each axis.
            const auto v = value();
            const auto colon = v ? v->find(':') : std::string_view::npos;
            const auto corner = [](std::string_view t) -> std::optional<Point> {
                Point q{};
                for (std::size_t k = 0; k < 3; ++k) {
                    const auto comma = t.find(',');
                    if ((k < 2) == (comma == std::string_view::npos)) return std::nullopt;
                    const auto x = number(t.substr(0, comma));
                    if (!x) return std::nullopt;
                    q[k] = *x * 1e-4;
                    if (k < 2) t = t.substr(comma + 1);
                }
                return q;
            };
            const auto a0 = colon == std::string_view::npos ? std::nullopt : corner(v->substr(0, colon));
            const auto a1 = colon == std::string_view::npos ? std::nullopt : corner(v->substr(colon + 1));
            if (!a0 || !a1 || !((*a0)[0] < (*a1)[0] && (*a0)[1] < (*a1)[1] && (*a0)[2] < (*a1)[2])) {
                return usage("--clip needs <x0>,<y0>,<z0>:<x1>,<y1>,<z1> (um) with each lo below its hi");
            }
            o.clip = Box{*a0, *a1};
        }
        else if (a == "--arrows") {
            const auto v = value();
            const auto n = v ? number(*v) : std::nullopt;
            if (!n || *n < 1 || *n > 100000 || *n != std::floor(*n)) return usage("--arrows needs a count from 1 to 100000");
            o.arrows = static_cast<std::size_t>(*n);
        }
        else if (a == "--from") { if (!flag_value(o.from)) return usage("--from needs a contact"); }
        else if (a == "--seeds") { if (!flag_value(o.seeds)) return usage("--seeds needs a segment or a plane"); }
        else if (a == "--cut") { if (!flag_value(o.cut)) return usage("--cut needs a segment"); }
        else if (a == "--compare") { if (!flag_value(o.compare)) return usage("--compare needs a deck or a result file"); }
        else if (a == "--lines") {
            const auto v = value();
            const auto n = v ? number(*v) : std::nullopt;
            if (!n || *n < 1 || *n > 10000 || *n != std::floor(*n)) return usage("--lines needs a count from 1 to 10000");
            o.lines = static_cast<std::size_t>(*n);
        }
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
        } else if (a == "--levels") {
            const auto v = value();
            if (!v) return usage("--levels needs <v>[,<v>...]");
            for (std::string_view rest = *v;;) {
                const auto comma = rest.find(',');
                const auto level = number(rest.substr(0, comma));
                if (!level) return usage("--levels needs <v>[,<v>...]");
                o.levels.push_back(*level);
                if (comma == std::string_view::npos) break;
                rest = rest.substr(comma + 1);
            }
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
            o.size_given = true;
            o.height = static_cast<int>(*h);
        } else {
            return usage("unknown option " + std::string(a));
        }
    }
    if (o.out.empty()) return usage("--out is required");
    if (!o.streamlines.empty()) {
        if (o.streamlines != "E" && o.streamlines != "J" && o.streamlines != "Jn" && o.streamlines != "Jp") {
            return usage("--streamlines is E, J, Jn or Jp");
        }
        if (o.field.empty()) o.field = o.streamlines;  // drawn over their own magnitude
    }
    if (!o.vectors.empty()) {
        if (o.vectors != "E" && o.vectors != "J" && o.vectors != "Jn" && o.vectors != "Jp") {
            return usage("--vectors is E, J, Jn or Jp");
        }
        if (o.field.empty()) o.field = o.vectors;
    }
    if (o.view.empty()) o.view = !o.slices.empty() ? "slice" : o.field.empty() ? "solid" : "field";
    const bool field_given = !o.field.empty();  // the exploded view draws regions without one
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
        // The doping and the mesh quality need no results (a field, profile or difference of the doping).
        const bool doping = o.field == "doping" || o.field == "ND" || o.field == "NA";
        const bool needs_results = o.view != "solid" && o.view != "mesh" && o.view != "quality" &&
                                   o.view != "quality-histogram" &&
                                   !((o.view == "field" || o.view == "profile" || o.view == "difference") && doping);
        if (needs_results) {
            if (auto r = load_outputs(*plan, l); !r) return fail(r.error());
        }
    }
    const mesh::Mesh& m = device->mesh();
    // Unit 45 (RQ4): a result that names its mesh must be of the device's; older results (digest 0)
    // fall back to the node count below.
    {
        const std::uint64_t digest = result_digest(l);
        if (digest != 0 && digest != mesh::graph_digest(m)) {
            return fail("the result was computed on another mesh: its mesh digest is not the device's");
        }
    }
    // Run B (Unit 26b, XD6): a result file read with run A's device, or a deck with its device and results.
    // A difference needs one mesh (the digest; for an older file, the node count).
    Loaded lb;
    std::optional<deck::Plan> plan_b;
    const device::Device* device_b = nullptr;
    bool same_mesh_b = false;
    if (!o.compare.empty()) {
        const std::string bext = o.compare.extension().string();
        if (bext == ".nitsweep" || bext == ".nittran" || bext == ".nitac") {
            if (auto r = load_result(o.compare, lb); !r) return fail(r.error());
            device_b = device;
            const std::uint64_t digest = result_digest(lb);
            const auto pb = points_of(lb);
            same_mesh_b = digest != 0 ? digest == mesh::graph_digest(m)
                                      : !pb.empty() && pb.front().fields.potential_V.size() == m.node_count();
        } else {
            const auto deck_b = deck::read(o.compare);
            if (!deck_b) return fail(deck_b.error().message);
            auto built = deck::build(*deck_b);
            if (!built) return fail(built.error().message);
            plan_b = std::move(*built);
            device_b = &plan_b->device;
            const bool doping = o.field == "doping" || o.field == "ND" || o.field == "NA";
            if (!(o.view == "difference" && doping)) {
                if (auto r = load_outputs(*plan_b, lb); !r) return fail(r.error());
            }
            same_mesh_b = mesh::graph_digest(device_b->mesh()) == mesh::graph_digest(m);
        }
    }
    const std::vector<results::BiasPoint> points_b = points_of(lb);
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
    if (o.clip) {
        if (m.dimension() != 3) return fail("a clipping box needs a 3D device");
        if (o.view != "solid" && o.view != "field") return fail("a clipping box applies to the solid and field views");
        view.clip = o.clip;
    }
    // A cut line (Unit 26b, XD5): plotted along by the profile and bands views, drawn on 2D and 3D device views.
    std::optional<CutLine> cut;
    if (!o.cut.empty()) {
        const bool device_view = o.view == "solid" || o.view == "mesh" || o.view == "field" || o.view == "slice" ||
                                 o.view == "iso" || o.view == "quality" || o.view == "difference";
        if (o.view != "profile" && o.view != "bands" && !(device_view && m.dimension() >= 2)) {
            return fail("a cut line applies to the profile and bands views and to 2D and 3D device views");
        }
        const auto ab = segment_um(o.cut, m.dimension());
        if (!ab) return usage("--cut needs <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] in the device's dimension");
        auto c = cut_line(*geometry, ab->first, ab->second);
        if (!c) return fail(c.error().message);
        cut = std::move(*c);
    }
    view.transparent_insulators = !o.opaque;
    view.projection = o.perspective ? Projection::perspective : Projection::orthographic;
    if (o.azimuth) view.azimuth_deg = *o.azimuth;
    if (o.elevation) view.elevation_deg = *o.elevation;
    std::optional<FieldScene> field_scene;
    std::optional<PlotScene> plot_scene;
    bool colored_arrows = false;  // 3D arrows on the solid, colored by their own magnitude
    std::vector<FieldScene> film;  // a filmstrip's frames
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
        auto v = field_values(*device, *geometry, **p, *field);
        if (!v) return std::unexpected(v.error().message);
        // A vector magnitude on a log scale (Unit 25b): without --range, the ten decades below its
        // largest value, everything under them (zero included) at the bottom color. Below that the
        // reconstructed currents are rounding: Jn at equilibrium spans 1e-42 to 1e-10 A/cm^2.
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
    const auto field_options = [&](const std::vector<double>& values) -> std::expected<FieldOptions, std::string> {
        FieldOptions f;
        f.view = view;
        f.map = *colormap;
        f.contours = o.contours;
        f.title = std::string(field_info(*field).title);
        // A magnitude that is zero everywhere (a quasi-static result's currents) is drawn on 0..1.
        const bool all_zero = is_vector_magnitude(*field) &&
                              std::ranges::none_of(values, [](double x) { return std::isfinite(x) && x != 0.0; });
        auto s = o.range    ? ColorScale::create(*scale_kind, o.range->first, o.range->second, 1e10)
                 : all_zero ? ColorScale::create(ScaleKind::linear, 0.0, 1.0, 1e10)
                            : ColorScale::fit(*scale_kind, values, 1e10);
        if (!s) return std::unexpected(s.error().message);
        f.scale = *s;
        return f;
    };

    if (o.view == "solid" || o.view == "mesh") {
        auto s = o.view == "solid" ? solid_view(*device, *geometry, view) : mesh_view(*device, *geometry, view);
        if (!s) return fail(s.error().message);
        field_scene = std::move(*s);
    } else if (!o.vectors.empty() && o.streamlines.empty() && m.dimension() == 3 && o.view == "field" &&
               o.field == o.vectors) {
        // 3D arrows over their own magnitude: on the solid device, colored by it.
        auto v = values_for();
        if (!v) return fail(v.error());
        auto f = field_options(*v);
        if (!f) return fail(f.error());
        auto s = solid_view(*device, *geometry, view);
        if (!s) return fail(s.error().message);
        s->title = f->title;
        s->colorbar = Colorbar{f->map, *f->scale, f->title};
        field_scene = std::move(*s);
        colored_arrows = true;
    } else if (!o.streamlines.empty() && m.dimension() == 3 && o.view == "field" && o.field == o.streamlines) {
        // 3D streamlines over their own magnitude: the lines, colored by it, in the device's outline.
        auto v = values_for();
        if (!v) return fail(v.error());
        auto f = field_options(*v);
        if (!f) return fail(f.error());
        auto s = outline_view(*device, *geometry, view);
        if (!s) return fail(s.error().message);
        s->title = f->title;
        s->colorbar = Colorbar{f->map, *f->scale, f->title};
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
    } else if (o.view == "iso") {
        if (o.levels.empty()) return usage("--view iso needs --levels");
        if (m.dimension() == 1) return fail("isosurfaces need a 2D or 3D device (a 1D field is a profile)");
        auto v = values_for();
        if (!v) return fail(v.error());
        auto f = field_options(*v);
        if (!f) return fail(f.error());
        auto s = iso_view(*device, *geometry, *v, o.levels, *f);
        if (!s) return fail(s.error().message);
        field_scene = std::move(*s);
    } else if (o.view == "filmstrip") {
        if (!l.transient) return fail("a filmstrip needs a transient run's results");
        if (m.dimension() == 1) return fail("a filmstrip needs a 2D or 3D device");
        if (points.empty()) return fail("the transient run stored no snapshots");
        // At most 12 frames: every k-th snapshot from the first.
        const std::size_t every = o.every > 0 ? o.every : (points.size() + 11) / 12;
        std::vector<FilmFrame> frames;
        for (std::size_t k = 0; k < points.size() && frames.size() < 12; k += every) {
            if (points[k].fields.potential_V.size() != m.node_count()) {
                return fail("the result's snapshots do not belong to the device's mesh");
            }
            auto v = field_values(*device, *geometry, points[k], *field);
            if (!v) return fail(v.error().message);
            frames.push_back({time_label(l.transient->snapshots[k].time_s), std::move(*v)});
        }
        // A vector magnitude on a log scale: the ten decades below the largest of all frames.
        if (is_vector_magnitude(*field) && *scale_kind == ScaleKind::log10 && !o.range) {
            double top = 0.0;
            for (const FilmFrame& f : frames) for (const double x : f.values) if (std::isfinite(x)) top = std::max(top, x);
            if (top > 0.0) {
                const double floor = std::pow(10.0, std::ceil(std::log10(top)) - 10.0);
                for (FilmFrame& f : frames) for (double& x : f.values) if (std::isfinite(x) && x < floor) x = floor;
            }
        }
        std::vector<double> all;
        for (const FilmFrame& f : frames) all.insert(all.end(), f.values.begin(), f.values.end());
        auto f = field_options(all);  // CQ2: one scale over every frame
        if (!f) return fail(f.error());
        auto s = filmstrip_scenes(*device, *geometry, frames, *scale_kind, *f);
        if (!s) return fail(s.error().message);
        film = std::move(*s);
    } else if (o.view == "exploded") {
        if (m.dimension() != 3) return fail("the exploded view needs a 3D device");
        const double gap = o.gap ? *o.gap * 1e-4 : default_gap(*device, o.axis);
        std::expected<FieldScene, base::Error> s;
        if (field_given) {
            auto v = values_for();
            if (!v) return fail(v.error());
            auto f = field_options(*v);
            if (!f) return fail(f.error());
            s = exploded_view(*device, *geometry, o.axis, gap, *v, *f);
        } else {
            FieldOptions f;
            f.view = view;
            s = exploded_view(*device, *geometry, o.axis, gap, {}, f);
        }
        if (!s) return fail(s.error().message);
        field_scene = std::move(*s);
    } else if (o.view == "profile" || o.view == "bands") {
        if (m.dimension() != 1 && !cut) return fail("--view " + o.view + " in 2D and 3D needs --cut");
        if (o.view == "bands") {
            auto p = pick();
            if (!p) return fail(p.error());
            auto s = cut ? cut_band_plot(m, *cut, (*p)->bands) : band_plot(m, (*p)->bands);
            if (!s) return fail(s.error().message);
            plot_scene = std::move(*s);
        } else {
            // A profile of the doping needs no results.
            auto v = values_for();
            if (!v) return fail(v.error());
            const AxisScale y = o.log ? AxisScale::log10 : AxisScale::linear;
            std::string title(field_info(*field).title);
            auto s = cut ? cut_plot(m, *cut, *v, std::move(title), y) : profile_plot(m, *v, std::move(title), y);
            if (!s) return fail(s.error().message);
            plot_scene = std::move(*s);
        }
    } else if (o.view == "iv" || o.view == "cv" || o.view == "current") {
        const auto contacts = device->contacts();
        std::vector<SeriesInput> series;
        std::string x_label, y_label;
        AxisScale x_scale = AxisScale::linear;
        // Run B's curves beside run A's (Unit 26b, XD6): the contacts chosen on A, found by name in B.
        const bool with_b = device_b != nullptr;
        const auto in_b = [&](const std::string& name) -> std::optional<std::size_t> {
            const std::size_t k = contact_index(*device_b, name, device_b->contacts().size());
            if (k == device_b->contacts().size()) return std::nullopt;
            return k;
        };
        const auto labelled = [&](std::string name, const char* run) { return with_b ? name + " (" + run + ")" : name; };
        if (o.view == "iv") {
            if (!l.sweep) return fail("an I-V plot needs a sweep's results");
            std::vector<std::vector<double>> biases;
            for (const auto& p : l.sweep->points) biases.push_back(p.bias_V);
            const std::size_t swept = contact_index(*device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(*device, o.sense, swept);
            auto c = analysis::current_curve(*l.sweep, swept, sense);
            if (!c) return fail(c.error().message);
            series.push_back({labelled("I(" + contacts[sense].name + ")", "A"), std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
            if (with_b) {
                if (!lb.sweep) return fail("run B holds no sweep to compare");
                const auto bs = in_b(contacts[swept].name), bn = in_b(contacts[sense].name);
                if (!bs) return fail("run B has no contact " + contacts[swept].name);
                if (!bn) return fail("run B has no contact " + contacts[sense].name);
                auto cb = analysis::current_curve(*lb.sweep, *bs, *bn);
                if (!cb) return fail(cb.error().message);
                series.push_back({labelled("I(" + contacts[sense].name + ")", "B"), std::move(*cb)});
            }
        } else if (o.view == "cv") {
            if (!l.small_signal) return fail("a C-V plot needs a small-signal run's results");
            std::vector<std::vector<double>> biases;
            for (const auto& p : l.small_signal->points) biases.push_back(p.dc.bias_V);
            const std::size_t swept = contact_index(*device, o.contact, swept_contact(biases));
            const std::size_t sense = contact_index(*device, o.sense, swept);
            const std::size_t f = o.frequency.value_or(0);
            auto c = analysis::capacitance_curve(*l.small_signal, f, sense, swept, swept);
            if (!c) return fail(c.error().message);
            const std::string name = "C at " + std::to_string(l.small_signal->frequency_Hz[f]) + " Hz";
            series.push_back({labelled(name, "A"), std::move(*c)});
            x_label = "V(" + contacts[swept].name + ") (V)";
            y_label = m.dimension() == 1 ? "C (F/cm\xC2\xB2)" : m.dimension() == 2 ? "C (F/cm)" : "C (F)";
            if (with_b) {
                if (!lb.small_signal || f >= lb.small_signal->frequency_Hz.size() ||
                    lb.small_signal->frequency_Hz[f] != l.small_signal->frequency_Hz[f]) {
                    return fail("run B holds no small-signal run at this frequency to compare");
                }
                const auto bs = in_b(contacts[swept].name), bn = in_b(contacts[sense].name);
                if (!bs) return fail("run B has no contact " + contacts[swept].name);
                if (!bn) return fail("run B has no contact " + contacts[sense].name);
                auto cb = analysis::capacitance_curve(*lb.small_signal, f, *bn, *bs, *bs);
                if (!cb) return fail(cb.error().message);
                series.push_back({labelled(name, "B"), std::move(*cb)});
            }
        } else {
            if (!l.transient) return fail("a current plot needs a transient run's results");
            if (with_b && !lb.transient) return fail("run B holds no transient run to compare");
            const auto curve_of = [](const results::Transient& t, std::size_t c) {
                std::vector<double> time, i;
                for (const auto& p : t.points) {
                    time.push_back(p.time_s);
                    i.push_back(p.terminal_current[c]);
                }
                return analysis::make_curve(std::move(time), std::move(i));
            };
            for (std::size_t k = 0; k < contacts.size(); ++k) {
                if (!o.contact.empty() && contacts[k].name != o.contact) continue;
                auto c = curve_of(*l.transient, k);
                if (!c) return fail(c.error().message);
                series.push_back({labelled(contacts[k].name, "A"), std::move(*c)});
                if (with_b) {
                    const auto bk = in_b(contacts[k].name);
                    if (!bk) return fail("run B has no contact " + contacts[k].name);
                    auto cb = curve_of(*lb.transient, *bk);
                    if (!cb) return fail(cb.error().message);
                    series.push_back({labelled(contacts[k].name, "B"), std::move(*cb)});
                }
            }
            x_label = "t (s)";
            y_label = "I (" + current_unit(m.dimension()) + ")";
        }
        auto s = curve_plot(series, x_label, y_label, x_scale, o.log ? AxisScale::log10 : AxisScale::linear);
        if (!s) return fail(s.error().message);
        plot_scene = std::move(*s);
    } else if (o.view == "quality" || o.view == "quality-histogram") {
        // Mesh quality (Unit 26b, XD7).
        auto q = element_quality(*geometry);
        if (!q) return fail(q.error().message);
        if (o.view == "quality-histogram") {
            auto s = quality_plot(*q, m.dimension());
            if (!s) return fail(s.error().message);
            plot_scene = std::move(*s);
        } else {
            FieldOptions f;
            f.view = view;
            f.map = *colormap;
            if (o.range) {
                auto s = ColorScale::create(ScaleKind::linear, o.range->first, o.range->second);
                if (!s) return fail(s.error().message);
                f.scale = *s;
            }
            auto s = quality_view(*device, *geometry, *q, o.slices, f);
            if (!s) return fail(s.error().message);
            field_scene = std::move(*s);
        }
    } else if (o.view == "difference") {
        // Run B - run A (Unit 26b, XD6), at the bias point of B equal to A's, on one mesh.
        if (!device_b) return fail("a difference map needs run B (--compare)");
        if (!same_mesh_b) return fail("run B was computed on another mesh: a difference needs one mesh");
        const bool doping = *field == Field::net_doping || *field == Field::donors || *field == Field::acceptors;
        const results::BiasPoint none;
        const results::BiasPoint* pa = &none;
        const results::BiasPoint* pb = &none;
        if (!doping) {
            auto p = pick();
            if (!p) return fail(p.error());
            pa = *p;
            const auto k = matching_point(*pa, points_b);
            if (!k) return fail("run B has no bias point at run A's biases");
            pb = &points_b[*k];
        }
        auto d = difference_field(*device, *device_b, *geometry, *pa, *pb, *field);
        if (!d) return fail(d.error().message);
        FieldOptions f;
        f.view = view;
        f.map = o.colormap == "turbo" ? Colormap::coolwarm : *colormap;  // XD6: cool-warm unless one is chosen
        f.contours = o.contours;
        f.title = d->title;
        if (o.range) {
            auto s = ColorScale::create(ScaleKind::linear, o.range->first, o.range->second);
            if (!s) return fail(s.error().message);
            f.scale = *s;
        } else {
            f.scale = difference_scale(d->values);
        }
        if (m.dimension() == 1) {
            auto p = profile_plot(m, d->values, d->title);
            if (!p) return fail(p.error().message);
            plot_scene = std::move(*p);
        } else {
            auto s = m.dimension() == 3 && !o.slices.empty() ? slice_view(*device, *geometry, d->values, o.slices, f)
                                                             : field_view(*device, *geometry, d->values, f);
            if (!s) return fail(s.error().message);
            field_scene = std::move(*s);
        }
    } else {
        return usage("unknown view " + o.view);
    }

    // Streamlines (Unit 25b, DS5-DS7) over the view.
    if (!o.streamlines.empty()) {
        if (!field_scene || m.dimension() < 2) {
            return fail("streamlines need a 2D or 3D device and a field, slice or iso view");
        }
        const VectorField vf = o.streamlines == "E"    ? VectorField::electric
                               : o.streamlines == "J"  ? VectorField::current
                               : o.streamlines == "Jn" ? VectorField::electron_current
                                                       : VectorField::hole_current;
        auto p = pick();
        if (!p) return fail(p.error());
        auto ev = element_vectors(*device, *geometry, **p, vf);
        if (!ev) return fail(ev.error().message);
        const std::vector<Point> nodal = streamline_field(*device, *geometry, *ev, vf);
        const Locator locator(*geometry);
        const auto contacts = device->contacts();
        std::vector<Point> seeds;
        if (!o.seeds.empty()) {
            const std::string_view s = o.seeds;
            if (s.size() > 2 && s[1] == '=' && std::string_view("xyz").find(s[0]) != std::string_view::npos) {
                const auto at = number(s.substr(2));
                if (!at) return usage("--seeds needs a segment or x=<um>, y=<um> or z=<um>");
                if (m.dimension() != 3) return usage("a seed plane needs a 3D device");
                Plane plane;
                plane.normal[static_cast<std::size_t>(s[0] - 'x')] = 1.0;
                plane.offset = *at * 1e-4;
                seeds = plane_seeds(*geometry, plane, o.lines);
                if (seeds.empty()) return fail("the seed plane misses the device");
            } else {
                const auto ab = segment_um(s, m.dimension());
                if (!ab) return usage("--seeds needs <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] in the device's dimension");
                seeds = segment_seeds(ab->first, ab->second, o.lines);
            }
        } else {
            // From a contact: --from, else for currents the one carrying the most, for E the first
            // gate or electrode, else the most biased.
            std::size_t c = 0;
            if (!o.from.empty()) {
                const device::Contact* named = device->find_contact(o.from);
                if (!named) return fail("there is no contact " + o.from);
                c = static_cast<std::size_t>(named - contacts.data());
            } else if (vf != VectorField::electric) {
                double best = -1.0;
                for (std::size_t k = 0; k < contacts.size(); ++k) {
                    double sum = 0.0;
                    for (const double i : contact_node_currents(*device, **p, k)) sum += i;
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
                vf == VectorField::electric ? std::vector<double>{} : contact_node_currents(*device, **p, c);
            seeds = contact_seeds(*device, locator, c, o.lines, currents);
        }
        std::vector<Streamline> lines;
        for (const Point& seed : seeds) lines.push_back(trace_streamline(locator, nodal, seed));
        if (m.dimension() == 3 && o.view == "field" && o.field == o.streamlines && field_scene->colorbar) {
            add_streamlines(*field_scene, lines, *field_scene->colorbar);
        } else {
            add_streamlines(*field_scene, lines, Rgba{0.05f, 0.05f, 0.05f, 0.85f});
        }
    }

    // Arrows (Unit 25c, DC1) over the view: 2D on a grid, 3D on the slices shown or the boundary.
    if (!o.vectors.empty()) {
        if (!field_scene || m.dimension() < 2) {
            return fail("arrows need a 2D or 3D device and a field, slice, iso or solid view");
        }
        const VectorField vf = o.vectors == "E"    ? VectorField::electric
                               : o.vectors == "J"  ? VectorField::current
                               : o.vectors == "Jn" ? VectorField::electron_current
                                                   : VectorField::hole_current;
        auto p = pick();
        if (!p) return fail(p.error());
        auto ev = element_vectors(*device, *geometry, **p, vf);
        if (!ev) return fail(ev.error().message);
        const std::vector<Point> nodal = streamline_field(*device, *geometry, *ev, vf);
        const Locator locator(*geometry);
        const std::size_t count = o.arrows > 0 ? o.arrows : (m.dimension() == 2 ? 300 : 400);
        // Currents: no arrows in the insulators, which carry none.
        std::vector<std::uint8_t> insulator;
        if (vf != VectorField::electric) {
            for (std::size_t i = 0; i < m.node_count(); ++i) insulator.push_back(device->is_insulator(static_cast<mesh::NodeId>(i)) ? 1 : 0);
        }
        double spacing = 0.0;
        std::vector<Arrow> arrows;
        if (m.dimension() == 2) {
            arrows = grid_arrows(locator, nodal, count, spacing, insulator);
        } else {
            std::vector<Polygon> surface;
            if (o.view == "slice") {
                for (const Plane& plane : o.slices) {
                    for (Polygon& q : geometry->slice(plane)) surface.push_back(std::move(q));
                }
            } else {
                surface = exposed_boundary(*device, *geometry);  // the contacts' slabs cover theirs
            }
            arrows = surface_arrows(locator, nodal, surface, count, spacing, insulator);
        }
        if (colored_arrows && field_scene->colorbar) {
            add_arrows(*field_scene, arrows, 0.8 * spacing, *field_scene->colorbar);
        } else {
            add_arrows(*field_scene, arrows, 0.8 * spacing, Rgba{0.05f, 0.05f, 0.05f, 1.0f});
        }
    }

    if (cut && field_scene) add_cut_line(*field_scene, *cut);

    // Draw and write.
    const Adapter adapter = o.adapter == "gpu" ? Adapter::hardware : o.adapter == "warp" ? Adapter::warp : Adapter::automatic;
    if (o.adapter != "gpu" && o.adapter != "warp" && o.adapter != "auto") return usage("unknown adapter " + o.adapter);
    auto renderer = Renderer::create({.adapter = adapter});
    if (!renderer) return fail(renderer.error().message);
    if (!film.empty()) {
        // A filmstrip: each frame drawn at the frame size (--size, else 640x480), tiled in rows of
        // --columns (4); under --frames each is also written alone.
        const int columns = static_cast<int>(o.columns > 0 ? std::min(o.columns, film.size()) : std::min<std::size_t>(4, film.size()));
        const int tw = o.size_given ? o.width : 640, th = o.size_given ? o.height : 480;
        if (!o.frames.empty()) {
            std::error_code ec;
            fs::create_directories(o.frames, ec);
        }
        std::vector<Image> tiles;
        for (std::size_t k = 0; k < film.size(); ++k) {
            auto t = renderer->draw(film[k], tw, th);
            if (!t) return fail(t.error().message);
            if (!o.frames.empty()) {
                char name[32];
                std::snprintf(name, sizeof name, "frame_%03zu.png", k);
                if (auto saved = save_png(*t, o.frames / name); !saved) return fail(saved.error().message);
            }
            tiles.push_back(std::move(*t));
        }
        auto strip = tile_images(tiles, columns);
        if (!strip) return fail(strip.error().message);
        if (auto saved = save_png(*strip, o.out); !saved) return fail(saved.error().message);
        std::printf("wrote %s (%dx%d, %zu frames, %s)\n", o.out.string().c_str(), strip->width, strip->height,
                    film.size(), renderer->adapter_name().c_str());
        return 0;
    }
    auto image = field_scene ? renderer->draw(*field_scene, o.width, o.height) : renderer->draw(*plot_scene, o.width, o.height);
    if (!image) return fail(image.error().message);
    if (auto saved = save_png(*image, o.out); !saved) return fail(saved.error().message);
    std::printf("wrote %s (%dx%d, %s)\n", o.out.string().c_str(), image->width, image->height,
                renderer->adapter_name().c_str());
    return 0;
}

const char* plot_usage_text() { return plot_usage; }

}  // namespace NiTCAD::cli
