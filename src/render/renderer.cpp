// The off-screen renderer (image.hpp): Direct3D 12 for triangles and lines, Direct2D through
// D3D11On12 for text, colorbars, axes and plots.
#include "NiTCAD/render/image.hpp"

#include "NiTCAD/render/axes3d.hpp"
#include "NiTCAD/render/color.hpp"

#include <windows.h>
#include <d2d1_3.h>
#include <d3d11on12.h>
#include <d3d12.h>
#include <dwrite_3.h>
#include <dxgi1_6.h>
#include <wrl/client.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <format>
#include <functional>
#include <limits>
#include <map>
#include <new>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "background_ps.h"
#include "composite_ps.h"
#include "fullscreen_vs.h"
#include "line_ps.h"
#include "line_vs.h"
#include "mesh_vs.h"
#include "opaque_ps.h"
#include "resolve_ps.h"
#include "shadow_ps.h"
#include "shadow_vs.h"
#include "transparent_ps.h"

namespace NiTCAD::render {

namespace {

using Microsoft::WRL::ComPtr;

constexpr int supersample = 2;
constexpr DXGI_FORMAT color_format = DXGI_FORMAT_B8G8R8A8_UNORM;
constexpr DXGI_FORMAT accum_format = DXGI_FORMAT_R16G16B16A16_FLOAT;
constexpr DXGI_FORMAT reveal_format = DXGI_FORMAT_R16_FLOAT;
constexpr DXGI_FORMAT depth_format = DXGI_FORMAT_D32_FLOAT;

base::Error device_error(std::string_view what, HRESULT hr) {
    return {base::ErrorCode::resource_exhausted,
            std::format("{} failed (HRESULT 0x{:08X})", what, static_cast<unsigned>(hr)),
            base::ErrorContext{.index = std::nullopt, .value = static_cast<double>(static_cast<unsigned>(hr))}};
}

#define NITCAD_HR(call, what)                                         \
    do {                                                              \
        const HRESULT hr_ = (call);                                   \
        if (FAILED(hr_)) return std::unexpected(device_error(what, hr_)); \
    } while (false)

std::string narrow(const wchar_t* w) {
    const int n = WideCharToMultiByte(CP_UTF8, 0, w, -1, nullptr, 0, nullptr, nullptr);
    std::string s(static_cast<std::size_t>(std::max(n, 1)) - 1, '\0');
    if (n > 1) WideCharToMultiByte(CP_UTF8, 0, w, -1, s.data(), n, nullptr, nullptr);
    return s;
}

std::wstring widen(std::string_view s) {
    if (s.empty()) return {};
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
    std::wstring w(static_cast<std::size_t>(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), w.data(), n);
    return w;
}

// --- camera -------------------------------------------------------------------------------------

using Mat4 = std::array<float, 16>;  // row-major, row vectors: p' = p M

struct View {
    Mat4 view_projection{};
    std::array<double, 3> right{}, up{}, forward{};
    std::array<double, 3> eye_position{};  // perspective only, um
    bool perspective = false;
};

double dotd(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
std::array<double, 3> crossd(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
std::array<double, 3> unitd(std::array<double, 3> a) {
    const double n = std::sqrt(dotd(a, a));
    return n > 0.0 ? std::array<double, 3>{a[0] / n, a[1] / n, a[2] / n} : std::array<double, 3>{0, 0, 1};
}

// The camera of a scene fitted to its bounds at this aspect (width / height) of the view.
View make_view(const FieldScene& s, double aspect) {
    const double az = s.camera.azimuth_deg * std::numbers::pi / 180.0;
    const double el = s.camera.elevation_deg * std::numbers::pi / 180.0;
    // Eye direction from the target: +y up; azimuth from +z toward +x.
    const std::array<double, 3> eye{std::cos(el) * std::sin(az), std::sin(el), std::cos(el) * std::cos(az)};
    View v;
    v.forward = {-eye[0], -eye[1], -eye[2]};
    std::array<double, 3> world_up{0.0, 1.0, 0.0};
    if (std::abs(dotd(world_up, v.forward)) > 0.999) world_up = {0.0, 0.0, -1.0};
    v.right = unitd(crossd(v.forward, world_up));  // +x to the right seen from +z
    v.up = crossd(v.right, v.forward);
    // Fit: the bounds' corners in view coordinates.
    const std::array<double, 3> lo{s.lo[0], s.lo[1], s.lo[2]}, hi{s.hi[0], s.hi[1], s.hi[2]};
    const std::array<double, 3> c{0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300, radius = 0.0;
    for (int k = 0; k < 8; ++k) {
        const std::array<double, 3> p{(k & 1) ? hi[0] : lo[0], (k & 2) ? hi[1] : lo[1], (k & 4) ? hi[2] : lo[2]};
        const std::array<double, 3> d{p[0] - c[0], p[1] - c[1], p[2] - c[2]};
        x0 = std::min(x0, dotd(d, v.right));
        x1 = std::max(x1, dotd(d, v.right));
        y0 = std::min(y0, dotd(d, v.up));
        y1 = std::max(y1, dotd(d, v.up));
        radius = std::max(radius, std::sqrt(dotd(d, d)));
    }
    radius = std::max(radius, 1e-6);
    // The window's zoom and pan (Unit 26): with zoom 1 and no pan these are the fitted values exactly
    // (a division by 1 and an addition of 0 change no bit).
    const double half_h =
        std::max({0.5 * (y1 - y0), 0.5 * (x1 - x0) / aspect, 1e-6}) / 0.92 / static_cast<double>(s.camera.zoom);
    const double half_w = half_h * aspect;
    const double cx = 0.5 * (x0 + x1) - static_cast<double>(s.camera.pan[0]) * half_w;
    const double cy = 0.5 * (y0 + y1) - static_cast<double>(s.camera.pan[1]) * half_h;
    // View coordinates: vx = right . (p - c) - cx, vy = up . (p - c) - cy, vz = forward . (p - c).
    const auto& r = v.right;
    const auto& u = v.up;
    const auto& f = v.forward;
    const double tx = -dotd(r, c) - cx, ty = -dotd(u, c) - cy, tz = -dotd(f, c);
    Mat4 m{};
    if (s.camera.projection == Projection::orthographic) {
        // ndc x = vx / half_w, y = vy / half_h, z = (vz + R) / (2R).
        const double R = radius * 1.05;
        for (int i = 0; i < 3; ++i) {
            m[static_cast<std::size_t>(i * 4 + 0)] = static_cast<float>(r[static_cast<std::size_t>(i)] / half_w);
            m[static_cast<std::size_t>(i * 4 + 1)] = static_cast<float>(u[static_cast<std::size_t>(i)] / half_h);
            m[static_cast<std::size_t>(i * 4 + 2)] = static_cast<float>(f[static_cast<std::size_t>(i)] / (2.0 * R));
        }
        m[12] = static_cast<float>(tx / half_w);
        m[13] = static_cast<float>(ty / half_h);
        m[14] = static_cast<float>((tz + R) / (2.0 * R));
        m[15] = 1.0f;
    } else {
        // The eye at distance D behind the fitted view plane; w = vz + D.
        const double t = std::tan(0.5 * s.camera.fov_deg * std::numbers::pi / 180.0);
        const double D = half_h / t + radius;
        const double n = std::max(D - radius * 1.05, 1e-3 * D), fz = D + radius * 1.05;
        const double sx = 1.0 / (t * aspect), sy = 1.0 / t;
        for (int i = 0; i < 3; ++i) {
            m[static_cast<std::size_t>(i * 4 + 0)] = static_cast<float>(r[static_cast<std::size_t>(i)] * sx);
            m[static_cast<std::size_t>(i * 4 + 1)] = static_cast<float>(u[static_cast<std::size_t>(i)] * sy);
            m[static_cast<std::size_t>(i * 4 + 2)] = static_cast<float>(f[static_cast<std::size_t>(i)] * fz / (fz - n));
            m[static_cast<std::size_t>(i * 4 + 3)] = static_cast<float>(f[static_cast<std::size_t>(i)]);
        }
        m[12] = static_cast<float>(tx * sx);
        m[13] = static_cast<float>(ty * sy);
        m[14] = static_cast<float>((tz + D) * fz / (fz - n) - fz * n / (fz - n));
        m[15] = static_cast<float>(tz + D);
        // The eye: the point with vx = vy = 0 and w = vz + D = 0.
        for (std::size_t i = 0; i < 3; ++i) v.eye_position[i] = c[i] + cx * r[i] + cy * u[i] - D * f[i];
        v.perspective = true;
    }
    v.view_projection = m;
    return v;
}

// The light from the upper left, over the viewer's shoulder (both looks).
std::array<double, 3> key_light(const View& v) {
    std::array<double, 3> L{};
    for (std::size_t i = 0; i < 3; ++i) L[i] = -v.forward[i] * 0.75 + v.up[i] * 0.55 - v.right[i] * 0.35;
    return unitd(L);
}

}  // namespace

FrameLighting frame_lighting(const FieldScene& s, double aspect) noexcept {
    FrameLighting out;
    if (s.look != Look::enhanced || s.dimension != 3) return out;
    const View v = make_view(s, aspect);
    const std::array<double, 3> L = key_light(v);
    // The fill from the lower right, on the other side of the eye from the key.
    std::array<double, 3> F{};
    for (std::size_t i = 0; i < 3; ++i) F[i] = -v.forward[i] * 0.6 - v.up[i] * 0.25 + v.right[i] * 0.75;
    F = unitd(F);
    for (std::size_t i = 0; i < 3; ++i) {
        out.key[i] = static_cast<float>(L[i]);
        out.fill[i] = static_cast<float>(F[i]);
        out.up[i] = static_cast<float>(v.up[i]);
        out.eye[i] = static_cast<float>(v.perspective ? v.eye_position[i] : -v.forward[i]);
    }
    out.perspective = v.perspective;
    out.enhanced = true;
    return out;
}

namespace {

// Pixel position of a world point in a view rectangle.
D2D1_POINT_2F project(const View& v, const Float3& p, const D2D1_RECT_F& rect) {
    const Mat4& m = v.view_projection;
    const float x = p[0] * m[0] + p[1] * m[4] + p[2] * m[8] + m[12];
    const float y = p[0] * m[1] + p[1] * m[5] + p[2] * m[9] + m[13];
    const float w = p[0] * m[3] + p[1] * m[7] + p[2] * m[11] + m[15];
    const float nx = x / w, ny = y / w;
    return {rect.left + (nx * 0.5f + 0.5f) * (rect.right - rect.left),
            rect.top + (0.5f - ny * 0.5f) * (rect.bottom - rect.top)};
}

// --- GPU data ----------------------------------------------------------------------------------

// Root CBVs need 256-byte alignment: the block is padded to a multiple of it.
struct FrameConstants {
    Mat4 view_projection;
    float light[3];
    float ambient;
    float viewport[2];
    std::uint32_t flags;
    float line_width;
    float colormap[256][4];
    float no_data[4];
    // The enhanced look (Unit 26f, KD2, KD4); zero for classic (see scene.hlsl).
    float eye[3];  // toward the viewer (orthographic) or the camera's position (perspective), world
    float style;   // 1 enhanced
    float up[3];
    float perspective;  // 1 when eye is a position
    float fill[3];
    float spare;
    float bg_top[4];
    float bg_bottom[4];
    float view_rect[4];
    float half_vector[3];
    float spare2;
    float shadow_box[4];    // Unit 26g (MD5): footprint x lo, z lo, x hi, z hi, um
    float shadow_plane[4];  // floor y (um), margin (um), alpha at the centre, unused
};
static_assert(sizeof(FrameConstants) % 256 == 0);
static_assert(sizeof(FrameConstants) == 4352);  // the HLSL block's size; the fields match scene.hlsl

struct MeshVertex {
    float position[3];
    float normal[3];
    float value;
    float color[4];
};

struct LineVertex {
    float position[3];
    float other[3];
    float side;
    float color[4];
};

D3D12_RESOURCE_BARRIER transition(ID3D12Resource* r, D3D12_RESOURCE_STATES from, D3D12_RESOURCE_STATES to) {
    D3D12_RESOURCE_BARRIER b{};
    b.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
    b.Transition.pResource = r;
    b.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
    b.Transition.StateBefore = from;
    b.Transition.StateAfter = to;
    return b;
}

// The rectangle of a field image the 3D or 2D view fills: below the title band, left of the colorbar
// panel, inside the 2D axes' margins (the whole image when that leaves too little).
D2D1_RECT_F field_view_rect(const FieldScene& s, int w, int h) {
    const float title_h = s.title.empty() ? 12.0f : 36.0f;
    const float bar_w = s.colorbar ? 130.0f : 12.0f;
    const bool axes = s.dimension == 2;
    D2D1_RECT_F view{axes ? 66.0f : 12.0f, title_h, static_cast<float>(w) - bar_w, static_cast<float>(h) - (axes ? 48.0f : 12.0f)};
    if (view.right - view.left < 16.0f || view.bottom - view.top < 16.0f) {
        view = {0.0f, 0.0f, static_cast<float>(w), static_cast<float>(h)};
    }
    return view;
}

// A scene's triangle and line draws: where each batch lies in the vertex buffers, and how it is drawn.
struct MeshRange {
    std::size_t first, count;
    std::uint32_t flags;
    bool transparent;
};
struct LineRange {
    std::size_t first, count;
    float width;
    bool top;
    bool box;  // in the axis box's buffer (Unit 26i: rebuilt with the camera)
};

// Unit 26i (ND1): a FieldScene's vertex and line buffers on the GPU, kept by a window target so that a
// camera move draws them again. `scene` and the sizes identify what they were made from; `hash` is
// the geometry's (Debug builds only, ND3).
struct ResidentGeometry {
    ComPtr<ID3D12Resource> mesh, lines;
    std::size_t mesh_vertices = 0, line_vertices = 0;
    std::vector<MeshRange> ranges;
    std::vector<LineRange> line_ranges;  // the scene's own, without the axis box
    const FieldScene* scene = nullptr;
    std::vector<std::size_t> sizes;  // each triangle batch's vertices, then each line batch's points
    std::uint64_t hash = 0;
};

// Unit 26i: text layouts kept across frames, by format and text. A layout is immutable once made, so a kept one
// draws exactly as a new one; the overlays' labels (colour bar, axis box, view cube) repeat from frame to frame,
// and making their layouts was a measurable part of an enhanced frame (KG6). The formats are the renderer's
// own, made once with it.
using LayoutCache = std::map<std::pair<IDWriteTextFormat*, std::wstring>, ComPtr<IDWriteTextLayout>>;
constexpr std::size_t layout_cache_limit = 2048;
std::vector<std::size_t> batch_sizes(const FieldScene& s) {
    std::vector<std::size_t> out;
    out.reserve(s.triangles.size() + s.lines.size() + 1);
    out.push_back(static_cast<std::size_t>(s.dimension));
    for (const TriangleBatch& b : s.triangles) out.push_back(b.vertices.size());
    for (const LineBatch& l : s.lines) out.push_back(l.points.size());
    return out;
}

// FNV-1a over everything render_field turns into vertices and draw ranges (ND3).
std::uint64_t geometry_hash(const FieldScene& s) {
    std::uint64_t h = 14695981039346656037ull;
    const auto mix = [&](const void* p, std::size_t n) {
        const auto* b = static_cast<const unsigned char*>(p);
        for (std::size_t i = 0; i < n; ++i) h = (h ^ b[i]) * 1099511628211ull;
    };
    const auto f = [&](float v) { mix(&v, sizeof v); };
    for (const TriangleBatch& b : s.triangles) {
        const std::uint8_t k[4] = {static_cast<std::uint8_t>(b.shading), b.lit, b.transparent, static_cast<std::uint8_t>(b.finish)};
        mix(k, sizeof k);
        for (const SceneVertex& x : b.vertices) {
            for (float v : x.position) f(v);
            for (float v : x.normal) f(v);
            f(x.value);
            f(x.color.r), f(x.color.g), f(x.color.b), f(x.color.a);
        }
    }
    for (const LineBatch& l : s.lines) {
        f(l.color.r), f(l.color.g), f(l.color.b), f(l.color.a), f(l.width_px);
        const std::uint8_t top = l.on_top;
        mix(&top, 1);
        for (const Float3& p : l.points) for (float v : p) f(v);
    }
    return h;
}

}  // namespace

// --- a window's swap chain (Unit 26) ---------------------------------------------------------------

struct WindowTarget::Impl {
    HWND window = nullptr;
    ComPtr<IDXGISwapChain3> swap;
    std::array<ComPtr<ID3D12Resource>, 2> buffers;
    int width = 0, height = 0;
    ResidentGeometry resident;  // Unit 26i: the last FieldScene's buffers (released after flush)
    // Waits for the renderer's queue: a swap chain is released only after its last present has run.
    std::function<void()> flush;
    ~Impl() {
        if (flush) flush();
    }
};

// --- the renderer -------------------------------------------------------------------------------

struct Renderer::Impl {
    std::string adapter_name;
    bool warp = false;
    bool debug = false;

    ComPtr<IDXGIFactory6> factory;
    ComPtr<ID3D12Device> device;
    ComPtr<ID3D12CommandQueue> queue;
    ComPtr<ID3D12CommandAllocator> allocator;
    ComPtr<ID3D12GraphicsCommandList> list;
    ComPtr<ID3D12Fence> fence;
    UINT64 fence_value = 0;
    HANDLE event = nullptr;
    ComPtr<ID3D12InfoQueue> info;

    ComPtr<ID3D12RootSignature> root;
    ComPtr<ID3D12PipelineState> opaque, transparent, line_depth, line_top, composite, resolve, background, shadow;

    ComPtr<ID3D11Device> d11;
    ComPtr<ID3D11DeviceContext> d11_context;
    ComPtr<ID3D11On12Device> on12;
    ComPtr<ID2D1Factory3> d2d_factory;
    ComPtr<ID2D1Device2> d2d_device;
    ComPtr<ID2D1DeviceContext2> d2d;
    ComPtr<IDWriteFactory> dwrite;
    ComPtr<IDWriteTextFormat> text, small_text, title_text;

    // Per-size targets.
    int width = 0, height = 0;
    ComPtr<ID3D12Resource> ssaa, depth, accum, reveal, final_target, readback;
    ComPtr<ID3D12DescriptorHeap> rtv_heap, dsv_heap, srv_heap;
    UINT rtv_size = 0, srv_size = 0;
    int repeat_frames = 0;    // frame_time: passes to repeat
    double repeat_ms = 0.0;
    ComPtr<ID3D11Resource> wrapped;
    ComPtr<ID2D1Bitmap1> d2d_target;
    D3D12_PLACED_SUBRESOURCE_FOOTPRINT footprint{};
    D3D12_RESOURCE_STATES final_state = D3D12_RESOURCE_STATE_RENDER_TARGET;

    ~Impl() {
        if (queue && fence) (void)wait();
        if (event) CloseHandle(event);
    }

    std::expected<void, base::Error> wait() {
        ++fence_value;
        NITCAD_HR(queue->Signal(fence.Get(), fence_value), "ID3D12CommandQueue::Signal");
        if (fence->GetCompletedValue() < fence_value) {
            NITCAD_HR(fence->SetEventOnCompletion(fence_value, event), "SetEventOnCompletion");
            WaitForSingleObject(event, INFINITE);
        }
        return {};
    }

    std::expected<void, base::Error> execute() {
        NITCAD_HR(list->Close(), "ID3D12GraphicsCommandList::Close");
        ID3D12CommandList* lists[] = {list.Get()};
        queue->ExecuteCommandLists(1, lists);
        if (auto w = wait(); !w) return w;
        NITCAD_HR(allocator->Reset(), "ID3D12CommandAllocator::Reset");
        NITCAD_HR(list->Reset(allocator.Get(), nullptr), "ID3D12GraphicsCommandList::Reset");
        return {};
    }

    std::expected<void, base::Error> create(const RendererOptions& o);
    std::expected<void, base::Error> create_pipelines();
    std::expected<void, base::Error> ensure_targets(int w, int h);
    std::expected<ComPtr<ID3D12Resource>, base::Error> upload(const void* data, std::size_t bytes);
    std::expected<Image, base::Error> read_back();
    std::expected<void, base::Error> begin_2d();
    std::expected<void, base::Error> end_2d();
    // Unit 26i: with `keep`, the vertex buffers are kept there, and with `reuse` drawn again from there when
    // they were made from this scene (ND1).
    std::expected<void, base::Error> render_field(const FieldScene& s, int w, int h, ResidentGeometry* keep = nullptr,
                                                  bool reuse = false);
    bool reused = false;  // the last render_field drew kept buffers
    ComPtr<ID3D12Resource> box_upload;  // Unit 26i: the axis box vertices, rewritten each frame
    std::size_t box_upload_bytes = 0;
    LayoutCache layouts;  // Unit 26i: the overlays' text layouts, kept across frames
    std::expected<void, base::Error> render_plot(const PlotScene& s, int w, int h);
    std::expected<void, base::Error> show(WindowTarget::Impl& target);
};

std::expected<void, base::Error> Renderer::Impl::create(const RendererOptions& o) {
    debug = o.debug_layer;
    if (debug) {
        ComPtr<ID3D12Debug> d;
        NITCAD_HR(D3D12GetDebugInterface(IID_PPV_ARGS(&d)), "D3D12GetDebugInterface (is Graphics Tools installed?)");
        d->EnableDebugLayer();
    }
    NITCAD_HR(CreateDXGIFactory2(debug ? DXGI_CREATE_FACTORY_DEBUG : 0, IID_PPV_ARGS(&factory)), "CreateDXGIFactory2");
    ComPtr<IDXGIAdapter1> adapter;
    if (o.adapter != Adapter::warp) {
        for (UINT i = 0;; ++i) {
            ComPtr<IDXGIAdapter1> a;
            if (factory->EnumAdapterByGpuPreference(i, DXGI_GPU_PREFERENCE_HIGH_PERFORMANCE, IID_PPV_ARGS(&a)) ==
                DXGI_ERROR_NOT_FOUND) {
                break;
            }
            DXGI_ADAPTER_DESC1 desc{};
            a->GetDesc1(&desc);
            if (desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE) continue;
            // The Microsoft Basic Render Driver (WARP behind a display adapter, as on CI runners) is not
            // flagged software but is no hardware adapter.
            if (desc.VendorId == 0x1414 && desc.DeviceId == 0x8c) continue;
            if (SUCCEEDED(D3D12CreateDevice(a.Get(), D3D_FEATURE_LEVEL_11_0, __uuidof(ID3D12Device), nullptr))) {
                adapter = a;
                break;
            }
        }
        if (!adapter && o.adapter == Adapter::hardware) {
            return std::unexpected(device_error("finding a Direct3D 12 hardware adapter", E_FAIL));
        }
    }
    if (!adapter) {
        NITCAD_HR(factory->EnumWarpAdapter(IID_PPV_ARGS(&adapter)), "IDXGIFactory4::EnumWarpAdapter");
        warp = true;
    }
    DXGI_ADAPTER_DESC1 desc{};
    adapter->GetDesc1(&desc);
    adapter_name = narrow(desc.Description);
    NITCAD_HR(D3D12CreateDevice(adapter.Get(), D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&device)), "D3D12CreateDevice");
    if (debug) {
        if (SUCCEEDED(device.As(&info))) {
            info->SetMessageCountLimit(10000);
            // D3D11On12's CreateWrappedResource probes the Direct3D 12 resource it wraps for D3D11
            // sharing properties, and the layer reports the probe of a resource D3D12 created as
            // an error. It is the only message the pipeline raises (measured), and it is not ours to
            // fix: filtered by its exact ID, every other message still counts.
            D3D12_MESSAGE_ID deny[] = {D3D12_MESSAGE_ID_REFLECTSHAREDPROPERTIES_INVALIDOBJECT};
            D3D12_INFO_QUEUE_FILTER filter{};
            filter.DenyList.NumIDs = 1;
            filter.DenyList.pIDList = deny;
            NITCAD_HR(info->PushStorageFilter(&filter), "ID3D12InfoQueue::PushStorageFilter");
        }
    }
    D3D12_COMMAND_QUEUE_DESC qd{};
    qd.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
    NITCAD_HR(device->CreateCommandQueue(&qd, IID_PPV_ARGS(&queue)), "CreateCommandQueue");
    NITCAD_HR(device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&allocator)),
              "CreateCommandAllocator");
    NITCAD_HR(device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, allocator.Get(), nullptr,
                                        IID_PPV_ARGS(&list)),
              "CreateCommandList");
    NITCAD_HR(device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&fence)), "CreateFence");
    event = CreateEventW(nullptr, FALSE, FALSE, nullptr);
    if (!event) return std::unexpected(device_error("CreateEvent", HRESULT_FROM_WIN32(GetLastError())));

    if (auto p = create_pipelines(); !p) return p;

    // Direct2D on the same queue through D3D11On12 (VQ1).
    IUnknown* queues[] = {queue.Get()};
    NITCAD_HR(D3D11On12CreateDevice(device.Get(), D3D11_CREATE_DEVICE_BGRA_SUPPORT, nullptr, 0, queues, 1, 0,
                                    &d11, &d11_context, nullptr),
              "D3D11On12CreateDevice");
    NITCAD_HR(d11.As(&on12), "QueryInterface ID3D11On12Device");
    D2D1_FACTORY_OPTIONS fo{};
    NITCAD_HR(D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, __uuidof(ID2D1Factory3), &fo,
                                reinterpret_cast<void**>(d2d_factory.GetAddressOf())),
              "D2D1CreateFactory");
    ComPtr<IDXGIDevice> dxgi_device;
    NITCAD_HR(d11.As(&dxgi_device), "QueryInterface IDXGIDevice");
    NITCAD_HR(d2d_factory->CreateDevice(dxgi_device.Get(), &d2d_device), "ID2D1Factory3::CreateDevice");
    NITCAD_HR(d2d_device->CreateDeviceContext(D2D1_DEVICE_CONTEXT_OPTIONS_NONE, &d2d), "CreateDeviceContext");
    d2d->SetTextAntialiasMode(D2D1_TEXT_ANTIALIAS_MODE_GRAYSCALE);
    NITCAD_HR(DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(IDWriteFactory),
                                  reinterpret_cast<IUnknown**>(dwrite.GetAddressOf())),
              "DWriteCreateFactory");
    NITCAD_HR(dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_NORMAL,
                                       DWRITE_FONT_STRETCH_NORMAL, 13.0f, L"en-us", &text),
              "CreateTextFormat");
    NITCAD_HR(dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_NORMAL,
                                       DWRITE_FONT_STRETCH_NORMAL, 11.5f, L"en-us", &small_text),
              "CreateTextFormat");
    NITCAD_HR(dwrite->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_SEMI_BOLD, DWRITE_FONT_STYLE_NORMAL,
                                       DWRITE_FONT_STRETCH_NORMAL, 15.0f, L"en-us", &title_text),
              "CreateTextFormat");
    return {};
}

std::expected<void, base::Error> Renderer::Impl::create_pipelines() {
    // Root signature: b0 as a root CBV, t0-t1 in a table.
    D3D12_DESCRIPTOR_RANGE range{D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 2, 0, 0, D3D12_DESCRIPTOR_RANGE_OFFSET_APPEND};
    D3D12_ROOT_PARAMETER params[2]{};
    params[0].ParameterType = D3D12_ROOT_PARAMETER_TYPE_CBV;
    params[0].Descriptor = {0, 0};
    params[0].ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;
    params[1].ParameterType = D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE;
    params[1].DescriptorTable = {1, &range};
    params[1].ShaderVisibility = D3D12_SHADER_VISIBILITY_PIXEL;
    D3D12_ROOT_SIGNATURE_DESC rd{2, params, 0, nullptr,
                                 D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT};
    ComPtr<ID3DBlob> blob, error;
    NITCAD_HR(D3D12SerializeRootSignature(&rd, D3D_ROOT_SIGNATURE_VERSION_1, &blob, &error),
              "D3D12SerializeRootSignature");
    NITCAD_HR(device->CreateRootSignature(0, blob->GetBufferPointer(), blob->GetBufferSize(), IID_PPV_ARGS(&root)),
              "CreateRootSignature");

    const D3D12_INPUT_ELEMENT_DESC mesh_layout[] = {
        {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
        {"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
        {"VALUE", 0, DXGI_FORMAT_R32_FLOAT, 0, 24, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
        {"COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 28, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0}};
    const D3D12_INPUT_ELEMENT_DESC line_layout[] = {
        {"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
        {"OTHER", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
        {"SIDE", 0, DXGI_FORMAT_R32_FLOAT, 0, 24, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
        {"COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 28, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0}};

    const auto base = [&]() {
        D3D12_GRAPHICS_PIPELINE_STATE_DESC d{};
        d.pRootSignature = root.Get();
        d.RasterizerState.FillMode = D3D12_FILL_MODE_SOLID;
        d.RasterizerState.CullMode = D3D12_CULL_MODE_NONE;
        d.RasterizerState.DepthClipEnable = TRUE;
        d.BlendState.RenderTarget[0].RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;
        d.SampleMask = UINT_MAX;
        d.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
        d.NumRenderTargets = 1;
        d.RTVFormats[0] = color_format;
        d.DSVFormat = depth_format;
        d.SampleDesc.Count = 1;
        d.DepthStencilState.DepthEnable = TRUE;
        d.DepthStencilState.DepthWriteMask = D3D12_DEPTH_WRITE_MASK_ALL;
        d.DepthStencilState.DepthFunc = D3D12_COMPARISON_FUNC_LESS_EQUAL;
        return d;
    };
    {
        auto d = base();
        d.InputLayout = {mesh_layout, 4};
        d.VS = {mesh_vs_bytecode, sizeof(mesh_vs_bytecode)};
        d.PS = {opaque_ps_bytecode, sizeof(opaque_ps_bytecode)};
        NITCAD_HR(device->CreateGraphicsPipelineState(&d, IID_PPV_ARGS(&opaque)), "CreateGraphicsPipelineState opaque");
    }
    {
        auto d = base();
        d.InputLayout = {mesh_layout, 4};
        d.VS = {mesh_vs_bytecode, sizeof(mesh_vs_bytecode)};
        d.PS = {transparent_ps_bytecode, sizeof(transparent_ps_bytecode)};
        d.NumRenderTargets = 2;
        d.RTVFormats[0] = accum_format;
        d.RTVFormats[1] = reveal_format;
        d.DepthStencilState.DepthWriteMask = D3D12_DEPTH_WRITE_MASK_ZERO;
        d.BlendState.IndependentBlendEnable = TRUE;
        auto& a = d.BlendState.RenderTarget[0];
        a.BlendEnable = TRUE;
        a.SrcBlend = a.DestBlend = D3D12_BLEND_ONE;
        a.SrcBlendAlpha = a.DestBlendAlpha = D3D12_BLEND_ONE;
        a.BlendOp = a.BlendOpAlpha = D3D12_BLEND_OP_ADD;
        a.RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;
        auto& r = d.BlendState.RenderTarget[1];
        r.BlendEnable = TRUE;
        r.SrcBlend = r.SrcBlendAlpha = D3D12_BLEND_ZERO;
        r.DestBlend = D3D12_BLEND_INV_SRC_COLOR;  // revealage *= 1 - alpha (its one channel)
        r.DestBlendAlpha = D3D12_BLEND_INV_SRC_ALPHA;  // alpha factors may not be colors
        r.BlendOp = r.BlendOpAlpha = D3D12_BLEND_OP_ADD;
        r.RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;
        NITCAD_HR(device->CreateGraphicsPipelineState(&d, IID_PPV_ARGS(&transparent)),
                  "CreateGraphicsPipelineState transparent");
    }
    for (const bool top : {false, true}) {
        auto d = base();
        d.InputLayout = {line_layout, 4};
        d.VS = {line_vs_bytecode, sizeof(line_vs_bytecode)};
        d.PS = {line_ps_bytecode, sizeof(line_ps_bytecode)};
        auto& a = d.BlendState.RenderTarget[0];
        a.BlendEnable = TRUE;
        a.SrcBlend = D3D12_BLEND_SRC_ALPHA;
        a.DestBlend = D3D12_BLEND_INV_SRC_ALPHA;
        a.SrcBlendAlpha = D3D12_BLEND_ONE;
        a.DestBlendAlpha = D3D12_BLEND_INV_SRC_ALPHA;
        a.BlendOp = a.BlendOpAlpha = D3D12_BLEND_OP_ADD;
        d.DepthStencilState.DepthWriteMask = D3D12_DEPTH_WRITE_MASK_ZERO;
        if (top) d.DepthStencilState.DepthEnable = FALSE;
        NITCAD_HR(device->CreateGraphicsPipelineState(&d, IID_PPV_ARGS(top ? &line_top : &line_depth)),
                  "CreateGraphicsPipelineState line");
    }
    {
        auto d = base();
        d.InputLayout = {nullptr, 0};
        d.VS = {fullscreen_vs_bytecode, sizeof(fullscreen_vs_bytecode)};
        d.PS = {composite_ps_bytecode, sizeof(composite_ps_bytecode)};
        d.DepthStencilState.DepthEnable = FALSE;
        d.DepthStencilState.DepthWriteMask = D3D12_DEPTH_WRITE_MASK_ZERO;
        d.DSVFormat = DXGI_FORMAT_UNKNOWN;
        auto& a = d.BlendState.RenderTarget[0];
        a.BlendEnable = TRUE;
        a.SrcBlend = D3D12_BLEND_SRC_ALPHA;
        a.DestBlend = D3D12_BLEND_INV_SRC_ALPHA;
        a.SrcBlendAlpha = D3D12_BLEND_ONE;
        a.DestBlendAlpha = D3D12_BLEND_INV_SRC_ALPHA;
        a.BlendOp = a.BlendOpAlpha = D3D12_BLEND_OP_ADD;
        NITCAD_HR(device->CreateGraphicsPipelineState(&d, IID_PPV_ARGS(&composite)),
                  "CreateGraphicsPipelineState composite");
        d.PS = {resolve_ps_bytecode, sizeof(resolve_ps_bytecode)};
        d.BlendState.RenderTarget[0].BlendEnable = FALSE;
        NITCAD_HR(device->CreateGraphicsPipelineState(&d, IID_PPV_ARGS(&resolve)),
                  "CreateGraphicsPipelineState resolve");
        d.PS = {background_ps_bytecode, sizeof(background_ps_bytecode)};
        NITCAD_HR(device->CreateGraphicsPipelineState(&d, IID_PPV_ARGS(&background)),
                  "CreateGraphicsPipelineState background");
        // Unit 26g (MD5): the ground shadow, blended over the gradient (alpha stays 1).
        d.VS = {shadow_vs_bytecode, sizeof(shadow_vs_bytecode)};
        d.PS = {shadow_ps_bytecode, sizeof(shadow_ps_bytecode)};
        auto& b = d.BlendState.RenderTarget[0];
        b.BlendEnable = TRUE;
        b.SrcBlend = D3D12_BLEND_SRC_ALPHA;
        b.DestBlend = D3D12_BLEND_INV_SRC_ALPHA;
        b.SrcBlendAlpha = D3D12_BLEND_ONE;
        b.DestBlendAlpha = D3D12_BLEND_INV_SRC_ALPHA;
        b.BlendOp = b.BlendOpAlpha = D3D12_BLEND_OP_ADD;
        NITCAD_HR(device->CreateGraphicsPipelineState(&d, IID_PPV_ARGS(&shadow)),
                  "CreateGraphicsPipelineState shadow");
    }
    return {};
}

std::expected<void, base::Error> Renderer::Impl::ensure_targets(int w, int h) {
    if (w == width && h == height && final_target) return {};
    if (auto r = wait(); !r) return r;
    d2d->SetTarget(nullptr);
    d2d_target.Reset();
    wrapped.Reset();
    d11_context->Flush();
    width = w;
    height = h;
    const D3D12_HEAP_PROPERTIES gpu{D3D12_HEAP_TYPE_DEFAULT, D3D12_CPU_PAGE_PROPERTY_UNKNOWN, D3D12_MEMORY_POOL_UNKNOWN, 1, 1};
    const auto texture = [&](int tw, int th, DXGI_FORMAT f, D3D12_RESOURCE_FLAGS flags, D3D12_RESOURCE_STATES state,
                             const D3D12_CLEAR_VALUE* clear, ComPtr<ID3D12Resource>& out) -> std::expected<void, base::Error> {
        D3D12_RESOURCE_DESC d{};
        d.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
        d.Width = static_cast<UINT64>(tw);
        d.Height = static_cast<UINT>(th);
        d.DepthOrArraySize = 1;
        d.MipLevels = 1;
        d.Format = f;
        d.SampleDesc.Count = 1;
        d.Flags = flags;
        NITCAD_HR(device->CreateCommittedResource(&gpu, D3D12_HEAP_FLAG_NONE, &d, state, clear, IID_PPV_ARGS(&out)),
                  "CreateCommittedResource texture");
        return {};
    };
    const int sw = w * supersample, sh = h * supersample;
    D3D12_CLEAR_VALUE white{color_format, {1.0f, 1.0f, 1.0f, 1.0f}};
    D3D12_CLEAR_VALUE zero{accum_format, {0.0f, 0.0f, 0.0f, 0.0f}};
    D3D12_CLEAR_VALUE one{reveal_format, {1.0f, 1.0f, 1.0f, 1.0f}};
    D3D12_CLEAR_VALUE far_depth{};
    far_depth.Format = depth_format;
    far_depth.DepthStencil.Depth = 1.0f;
    if (auto r = texture(sw, sh, color_format, D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET,
                         D3D12_RESOURCE_STATE_RENDER_TARGET, &white, ssaa); !r) return r;
    if (auto r = texture(sw, sh, depth_format, D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL,
                         D3D12_RESOURCE_STATE_DEPTH_WRITE, &far_depth, depth); !r) return r;
    if (auto r = texture(sw, sh, accum_format, D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET,
                         D3D12_RESOURCE_STATE_RENDER_TARGET, &zero, accum); !r) return r;
    if (auto r = texture(sw, sh, reveal_format, D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET,
                         D3D12_RESOURCE_STATE_RENDER_TARGET, &one, reveal); !r) return r;
    if (auto r = texture(w, h, color_format, D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET,
                         D3D12_RESOURCE_STATE_RENDER_TARGET, &white, final_target); !r) return r;
    final_state = D3D12_RESOURCE_STATE_RENDER_TARGET;

    D3D12_RESOURCE_DESC fd = final_target->GetDesc();
    UINT64 total = 0;
    device->GetCopyableFootprints(&fd, 0, 1, 0, &footprint, nullptr, nullptr, &total);
    const D3D12_HEAP_PROPERTIES rb{D3D12_HEAP_TYPE_READBACK, D3D12_CPU_PAGE_PROPERTY_UNKNOWN, D3D12_MEMORY_POOL_UNKNOWN, 1, 1};
    D3D12_RESOURCE_DESC bd{};
    bd.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    bd.Width = total;
    bd.Height = 1;
    bd.DepthOrArraySize = 1;
    bd.MipLevels = 1;
    bd.SampleDesc.Count = 1;
    bd.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    NITCAD_HR(device->CreateCommittedResource(&rb, D3D12_HEAP_FLAG_NONE, &bd, D3D12_RESOURCE_STATE_COPY_DEST, nullptr,
                                              IID_PPV_ARGS(&readback)),
              "CreateCommittedResource readback");

    D3D12_DESCRIPTOR_HEAP_DESC hd{D3D12_DESCRIPTOR_HEAP_TYPE_RTV, 4, D3D12_DESCRIPTOR_HEAP_FLAG_NONE, 0};
    NITCAD_HR(device->CreateDescriptorHeap(&hd, IID_PPV_ARGS(&rtv_heap)), "CreateDescriptorHeap RTV");
    hd = {D3D12_DESCRIPTOR_HEAP_TYPE_DSV, 1, D3D12_DESCRIPTOR_HEAP_FLAG_NONE, 0};
    NITCAD_HR(device->CreateDescriptorHeap(&hd, IID_PPV_ARGS(&dsv_heap)), "CreateDescriptorHeap DSV");
    hd = {D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV, 3, D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE, 0};
    NITCAD_HR(device->CreateDescriptorHeap(&hd, IID_PPV_ARGS(&srv_heap)), "CreateDescriptorHeap SRV");
    rtv_size = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);
    srv_size = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
    auto rtv = rtv_heap->GetCPUDescriptorHandleForHeapStart();
    for (ID3D12Resource* r : {ssaa.Get(), accum.Get(), reveal.Get(), final_target.Get()}) {
        device->CreateRenderTargetView(r, nullptr, rtv);
        rtv.ptr += rtv_size;
    }
    device->CreateDepthStencilView(depth.Get(), nullptr, dsv_heap->GetCPUDescriptorHandleForHeapStart());
    auto srv = srv_heap->GetCPUDescriptorHandleForHeapStart();
    for (ID3D12Resource* r : {accum.Get(), reveal.Get(), ssaa.Get()}) {
        device->CreateShaderResourceView(r, nullptr, srv);
        srv.ptr += srv_size;
    }

    D3D11_RESOURCE_FLAGS f11{};
    f11.BindFlags = D3D11_BIND_RENDER_TARGET;
    NITCAD_HR(on12->CreateWrappedResource(final_target.Get(), &f11, D3D12_RESOURCE_STATE_RENDER_TARGET,
                                          D3D12_RESOURCE_STATE_COPY_SOURCE, IID_PPV_ARGS(&wrapped)),
              "ID3D11On12Device::CreateWrappedResource");
    ComPtr<IDXGISurface> surface;
    NITCAD_HR(wrapped.As(&surface), "QueryInterface IDXGISurface");
    const D2D1_BITMAP_PROPERTIES1 bp = D2D1::BitmapProperties1(
        D2D1_BITMAP_OPTIONS_TARGET | D2D1_BITMAP_OPTIONS_CANNOT_DRAW,
        D2D1::PixelFormat(color_format, D2D1_ALPHA_MODE_PREMULTIPLIED), 96.0f, 96.0f);
    NITCAD_HR(d2d->CreateBitmapFromDxgiSurface(surface.Get(), &bp, &d2d_target), "CreateBitmapFromDxgiSurface");
    return {};
}

std::expected<ComPtr<ID3D12Resource>, base::Error> Renderer::Impl::upload(const void* data, std::size_t bytes) {
    const D3D12_HEAP_PROPERTIES up{D3D12_HEAP_TYPE_UPLOAD, D3D12_CPU_PAGE_PROPERTY_UNKNOWN, D3D12_MEMORY_POOL_UNKNOWN, 1, 1};
    D3D12_RESOURCE_DESC d{};
    d.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    d.Width = std::max<std::size_t>(bytes, 256);
    d.Height = 1;
    d.DepthOrArraySize = 1;
    d.MipLevels = 1;
    d.SampleDesc.Count = 1;
    d.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    ComPtr<ID3D12Resource> r;
    NITCAD_HR(device->CreateCommittedResource(&up, D3D12_HEAP_FLAG_NONE, &d, D3D12_RESOURCE_STATE_GENERIC_READ,
                                              nullptr, IID_PPV_ARGS(&r)),
              "CreateCommittedResource upload");
    void* p = nullptr;
    const D3D12_RANGE none{0, 0};
    NITCAD_HR(r->Map(0, &none, &p), "ID3D12Resource::Map");
    if (bytes) std::memcpy(p, data, bytes);
    r->Unmap(0, nullptr);
    return r;
}

std::expected<void, base::Error> Renderer::Impl::begin_2d() {
    ID3D11Resource* res[] = {wrapped.Get()};
    on12->AcquireWrappedResources(res, 1);
    d2d->SetTarget(d2d_target.Get());
    d2d->BeginDraw();
    d2d->SetTransform(D2D1::Matrix3x2F::Identity());
    return {};
}

std::expected<void, base::Error> Renderer::Impl::end_2d() {
    const HRESULT hr = d2d->EndDraw();
    d2d->SetTarget(nullptr);
    ID3D11Resource* res[] = {wrapped.Get()};
    on12->ReleaseWrappedResources(res, 1);
    d11_context->Flush();
    if (FAILED(hr)) return std::unexpected(device_error("ID2D1DeviceContext::EndDraw", hr));
    final_state = D3D12_RESOURCE_STATE_COPY_SOURCE;
    return {};
}

std::expected<Image, base::Error> Renderer::Impl::read_back() {
    if (final_state != D3D12_RESOURCE_STATE_COPY_SOURCE) {
        const auto b = transition(final_target.Get(), final_state, D3D12_RESOURCE_STATE_COPY_SOURCE);
        list->ResourceBarrier(1, &b);
        final_state = D3D12_RESOURCE_STATE_COPY_SOURCE;
    }
    D3D12_TEXTURE_COPY_LOCATION dst{readback.Get(), D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT, {}};
    dst.PlacedFootprint = footprint;
    D3D12_TEXTURE_COPY_LOCATION src{final_target.Get(), D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX, {}};
    src.SubresourceIndex = 0;
    list->CopyTextureRegion(&dst, 0, 0, 0, &src, nullptr);
    const auto back = transition(final_target.Get(), D3D12_RESOURCE_STATE_COPY_SOURCE, D3D12_RESOURCE_STATE_RENDER_TARGET);
    list->ResourceBarrier(1, &back);
    final_state = D3D12_RESOURCE_STATE_RENDER_TARGET;
    if (auto r = execute(); !r) return std::unexpected(std::move(r.error()));
    Image img{width, height, std::vector<std::uint8_t>(static_cast<std::size_t>(width) * height * 4)};
    void* p = nullptr;

    NITCAD_HR(readback->Map(0, nullptr, &p), "ID3D12Resource::Map readback");  // the whole buffer
    const auto* bytes = static_cast<const std::uint8_t*>(p);
    for (int y = 0; y < height; ++y) {
        const std::uint8_t* row = bytes + static_cast<std::size_t>(footprint.Footprint.RowPitch) * y;
        std::uint8_t* out = img.rgba.data() + static_cast<std::size_t>(width) * 4 * y;
        for (int x = 0; x < width; ++x) {
            out[4 * x + 0] = row[4 * x + 2];
            out[4 * x + 1] = row[4 * x + 1];
            out[4 * x + 2] = row[4 * x + 0];
            out[4 * x + 3] = 255;
        }
    }
    const D3D12_RANGE written{0, 0};
    readback->Unmap(0, &written);
    return img;
}

namespace {

// --- Direct2D helpers ---------------------------------------------------------------------------

D2D1_COLOR_F d2(const Rgb& c, float a = 1.0f) { return {c.r, c.g, c.b, a}; }


struct Ink {
    ID2D1DeviceContext2* dc;
    IDWriteFactory* dw;
    LayoutCache* cache;
    ComPtr<ID2D1SolidColorBrush> brush;
    explicit Ink(ID2D1DeviceContext2* c, IDWriteFactory* w, LayoutCache* kept = nullptr) : dc(c), dw(w), cache(kept) {
        dc->CreateSolidColorBrush(D2D1::ColorF(0, 0, 0), &brush);
    }
    ID2D1SolidColorBrush* color(const D2D1_COLOR_F& c) {
        brush->SetColor(c);
        return brush.Get();
    }
    // The layout of the text in the format, kept in the cache when there is one; null on failure.
    ComPtr<IDWriteTextLayout> layout_of(std::string_view s, IDWriteTextFormat* f) {
        std::wstring w = widen(s);
        if (cache) {
            if (const auto it = cache->find({f, w}); it != cache->end()) return it->second;
        }
        ComPtr<IDWriteTextLayout> layout;
        if (FAILED(dw->CreateTextLayout(w.c_str(), static_cast<UINT32>(w.size()), f, 4096.0f, 1024.0f, &layout))) return nullptr;
        if (cache) {
            if (cache->size() >= layout_cache_limit) cache->clear();
            cache->emplace(std::pair{f, std::move(w)}, layout);
        }
        return layout;
    }
    // Text anchored at (x, y): align 0 left, 1 center, 2 right; valign 0 top, 1 middle, 2 bottom.
    void text(std::string_view s, IDWriteTextFormat* f, float x, float y, int align, int valign,
              const D2D1_COLOR_F& c = {0.12f, 0.14f, 0.16f, 1.0f}) {
        const ComPtr<IDWriteTextLayout> layout = layout_of(s, f);
        if (!layout) return;
        DWRITE_TEXT_METRICS m{};
        layout->GetMetrics(&m);
        const float left = x - (align == 1 ? 0.5f * m.width : align == 2 ? m.width : 0.0f);
        const float top = y - (valign == 1 ? 0.5f * m.height : valign == 2 ? m.height : 0.0f);
        dc->DrawTextLayout({left, top}, layout.Get(), color(c));
    }
    float width_of(std::string_view s, IDWriteTextFormat* f) {
        const ComPtr<IDWriteTextLayout> layout = layout_of(s, f);
        if (!layout) return 0.0f;
        DWRITE_TEXT_METRICS m{};
        layout->GetMetrics(&m);
        return m.width;
    }
};

// The enhanced 3D look (Unit 26g, MD7): a lighter frame, and ticks inside the bar's right edge as well.
void draw_colorbar(Ink& ink, IDWriteTextFormat* text, IDWriteTextFormat* minor_text, const Colorbar& bar,
                   const D2D1_RECT_F& area, bool enhanced) {
    const float bw = 16.0f;
    const float top = area.top + 28.0f, bottom = area.bottom - 12.0f;
    const float left = area.left + 6.0f;
    constexpr int steps = 192;
    for (int k = 0; k < steps; ++k) {
        const double t = (k + 0.5) / steps;
        const float y1 = bottom - (bottom - top) * static_cast<float>(k) / steps;
        const float y0 = bottom - (bottom - top) * static_cast<float>(k + 1) / steps;
        ink.dc->FillRectangle({left, y0 - 0.5f, left + bw, y1 + 0.5f}, ink.color(d2(sample(bar.map, t))));
    }
    const D2D1_COLOR_F frame = enhanced ? D2D1_COLOR_F{0.45f, 0.47f, 0.5f, 1.0f} : D2D1_COLOR_F{0.3f, 0.32f, 0.35f, 1.0f};
    ink.dc->DrawRectangle({left, top, left + bw, bottom}, ink.color(frame), 1.0f);
    for (const Tick& t : bar.scale.ticks()) {
        const float y = bottom - (bottom - top) * static_cast<float>(t.position);
        const float len = t.major ? 5.0f : 3.0f;
        ink.dc->DrawLine({left + bw, y}, {left + bw + len, y}, ink.color(frame), 1.0f);
        if (enhanced) ink.dc->DrawLine({left + bw - len, y}, {left + bw, y}, ink.color(frame), 1.0f);
        if (t.major && !t.label.empty()) ink.text(t.label, minor_text, left + bw + 8.0f, y, 0, 1);
    }
    ink.text(bar.title, text, left, area.top + 6.0f, 0, 0);
}

}  // namespace

std::expected<void, base::Error> Renderer::Impl::render_field(const FieldScene& s, int w, int h, ResidentGeometry* keep,
                                                              bool reuse) {
    if (auto r = ensure_targets(w, h); !r) return std::unexpected(std::move(r.error()));
    // Layout: title band, colorbar panel on the right, 2D axes margins.
    const float title_h = s.title.empty() ? 12.0f : 36.0f;
    const float bar_w = s.colorbar ? 130.0f : 12.0f;
    const bool axes = s.dimension == 2;
    const D2D1_RECT_F view = field_view_rect(s, w, h);
    const double aspect = (view.right - view.left) / (view.bottom - view.top);
    const View v = make_view(s, aspect);

    // Vertex data: kept by a window target (Unit 26i, ND1), else this frame's own.
    ResidentGeometry own;
    ResidentGeometry& g = keep ? *keep : own;
    reused = reuse && keep && g.mesh && g.lines && g.scene == &s && g.sizes == batch_sizes(s);
#ifndef NDEBUG
    // ND3: the caller promised an unchanged geometry; a missed full present would draw a stale one.
    if (reused && g.hash != geometry_hash(s)) {
        std::fputs("NiTCAD render: present_camera on a scene whose geometry changed (Unit 26i, ND3)\n", stderr);
        std::abort();
    }
#endif
    // The staging buffers live until the frame has executed.
    ComPtr<ID3D12Resource> mesh_staging, line_staging, box_buffer;
    const auto resident = [&](const ComPtr<ID3D12Resource>& staging, std::size_t bytes)
        -> std::expected<ComPtr<ID3D12Resource>, base::Error> {
        // The vertices live in default-heap buffers: one copy, then read at GPU speed (a window
        // target keeps them across frames, Unit 26i).
        const D3D12_HEAP_PROPERTIES gpu{D3D12_HEAP_TYPE_DEFAULT, D3D12_CPU_PAGE_PROPERTY_UNKNOWN, D3D12_MEMORY_POOL_UNKNOWN, 1, 1};
        D3D12_RESOURCE_DESC d = staging->GetDesc();
        ComPtr<ID3D12Resource> r;
        NITCAD_HR(device->CreateCommittedResource(&gpu, D3D12_HEAP_FLAG_NONE, &d, D3D12_RESOURCE_STATE_COPY_DEST, nullptr,
                                                  IID_PPV_ARGS(&r)),
                  "CreateCommittedResource vertices");
        if (bytes) list->CopyBufferRegion(r.Get(), 0, staging.Get(), 0, bytes);
        const auto b = transition(r.Get(), D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_VERTEX_AND_CONSTANT_BUFFER);
        list->ResourceBarrier(1, &b);
        return r;
    };
    if (!reused) {
        g = ResidentGeometry{};  // the old buffers go; the queue is idle (wait() before every frame)
        std::vector<MeshVertex> mesh;
        for (const TriangleBatch& b : s.triangles) {
            const std::size_t first = mesh.size();
            for (const SceneVertex& x : b.vertices) {
                mesh.push_back({{x.position[0], x.position[1], x.position[2]},
                                {x.normal[0], x.normal[1], x.normal[2]},
                                x.value,
                                {x.color.r, x.color.g, x.color.b, x.color.a}});
            }
            // Bits 3-4, the finish (Unit 26g, MD1): 1 metal, 2 glass; read only by the enhanced look.
            const std::uint32_t flags = (b.lit && s.dimension == 3 ? 1u : 0u) | (b.shading == Shading::field ? 2u : 0u) |
                                        (static_cast<std::uint32_t>(b.finish) << 3);
            g.ranges.push_back({first, b.vertices.size(), flags, b.transparent});
        }
        std::vector<LineVertex> lines;
        for (const LineBatch& l : s.lines) {
            const std::size_t first = lines.size();
            for (std::size_t i = 0; i + 1 < l.points.size(); i += 2) {
                const Float3& a = l.points[i];
                const Float3& b = l.points[i + 1];
                const auto vert = [&](const Float3& p, const Float3& q, float side) {
                    lines.push_back({{p[0], p[1], p[2]}, {q[0], q[1], q[2]}, side, {l.color.r, l.color.g, l.color.b, l.color.a}});
                };
                // The quad a-, a+, b+, b- as two triangles. At b the direction to "other" is reversed,
                // so its geometric + side has side = -1.
                vert(a, b, -1.0f);
                vert(a, b, 1.0f);
                vert(b, a, -1.0f);
                vert(a, b, -1.0f);
                vert(b, a, -1.0f);
                vert(b, a, 1.0f);
            }
            g.line_ranges.push_back({first, lines.size() - first, l.width_px, l.on_top, false});
        }
        auto vb = upload(mesh.data(), mesh.size() * sizeof(MeshVertex));
        if (!vb) return std::unexpected(std::move(vb.error()));
        auto lb = upload(lines.data(), lines.size() * sizeof(LineVertex));
        if (!lb) return std::unexpected(std::move(lb.error()));
        mesh_staging = *vb;
        line_staging = *lb;
        auto vbr = resident(mesh_staging, mesh.size() * sizeof(MeshVertex));
        if (!vbr) return std::unexpected(std::move(vbr.error()));
        auto lbr = resident(line_staging, lines.size() * sizeof(LineVertex));
        if (!lbr) return std::unexpected(std::move(lbr.error()));
        g.mesh = std::move(*vbr);
        g.lines = std::move(*lbr);
        g.mesh_vertices = mesh.size();
        g.line_vertices = lines.size();
        if (keep) {
            g.scene = &s;
            g.sizes = batch_sizes(s);
#ifndef NDEBUG
            g.hash = geometry_hash(s);
#endif
        }
    }
    const std::vector<MeshRange>& ranges = g.ranges;
    std::vector<LineRange> line_ranges = g.line_ranges;
    // The enhanced look (Unit 26f): its lights, and the gradient behind a 3D view (KD2, KD4). Its
    // feature edges are already among the scene's lines (set_look).
    const FrameLighting lighting = frame_lighting(s, aspect);
    const bool enhanced = lighting.enhanced;
    // Unit 26g (MD5): not on clipped or exploded views, nor when looking up from below.
    const bool shadowed = enhanced && s.ground_shadow && s.camera.elevation_deg >= -5.0f;
    // Unit 26h (JD1): the axis box, its edges depth-tested with the scene's lines, its ticks on Direct2D.
    // It follows the camera, so it is rebuilt every frame, in a buffer of its own (Unit 26i).
    const AxisBox box = enhanced && s.show_axis_box ? axis_box(s, w, h) : AxisBox{};
    std::vector<LineVertex> box_lines;
    if (!box.edges.empty()) {
        const float c[4] = {0.55f, 0.58f, 0.62f, 1.0f};
        for (const auto& [a, b] : box.edges) {
            const auto vert = [&](const Float3& p, const Float3& q, float side) {
                box_lines.push_back({{p[0], p[1], p[2]}, {q[0], q[1], q[2]}, side, {c[0], c[1], c[2], c[3]}});
            };
            vert(a, b, -1.0f);
            vert(a, b, 1.0f);
            vert(b, a, -1.0f);
            vert(a, b, -1.0f);
            vert(b, a, -1.0f);
            vert(b, a, 1.0f);
        }
        line_ranges.push_back({0, box_lines.size(), 1.0f, false, true});
        // One upload-heap buffer kept by the renderer and rewritten each frame (the queue is idle between
        // frames): a fresh one per frame cost the enhanced draw some 5% (KG6).
        const std::size_t bytes = box_lines.size() * sizeof(LineVertex);
        if (!box_upload || box_upload_bytes < bytes) {
            auto bb = upload(box_lines.data(), bytes);  // at least 256 bytes (upload)
            if (!bb) return std::unexpected(std::move(bb.error()));
            box_upload = std::move(*bb);
            box_upload_bytes = std::max<std::size_t>(bytes, 256);
        }
        void* p = nullptr;
        const D3D12_RANGE none{0, 0};
        NITCAD_HR(box_upload->Map(0, &none, &p), "ID3D12Resource::Map axis box");
        std::memcpy(p, box_lines.data(), bytes);
        box_upload->Unmap(0, nullptr);
        box_buffer = box_upload;  // read from the upload heap: a few dozen vertices
    }

    // Constants: one block per draw.
    FrameConstants base{};
    base.view_projection = v.view_projection;
    {
        const std::array<double, 3> L = key_light(v);
        base.light[0] = static_cast<float>(L[0]);
        base.light[1] = static_cast<float>(L[1]);
        base.light[2] = static_cast<float>(L[2]);
    }
    base.ambient = 0.38f;
    base.viewport[0] = (view.right - view.left) * supersample;
    base.viewport[1] = (view.bottom - view.top) * supersample;
    if (enhanced) {
        for (std::size_t i = 0; i < 3; ++i) {
            base.eye[i] = lighting.eye[i];
            base.up[i] = lighting.up[i];
            base.fill[i] = lighting.fill[i];
        }
        base.perspective = lighting.perspective ? 1.0f : 0.0f;
        base.style = 1.0f;
        const float top[4] = {0.86f, 0.89f, 0.93f, 1.0f};
        const float bottom[4] = {0.99f, 0.99f, 1.0f, 1.0f};
        std::memcpy(base.bg_top, top, sizeof top);
        std::memcpy(base.bg_bottom, bottom, sizeof bottom);
        base.view_rect[0] = view.left * supersample;
        base.view_rect[1] = view.top * supersample;
        base.view_rect[2] = (view.right - view.left) * supersample;
        base.view_rect[3] = (view.bottom - view.top) * supersample;
        // The ground shadow (MD5): 4% of the height below the bounds, fading out one twelfth of the
        // largest extent outside the footprint.
        const float extent = std::max({s.hi[0] - s.lo[0], s.hi[1] - s.lo[1], s.hi[2] - s.lo[2]});
        base.shadow_box[0] = s.lo[0];
        base.shadow_box[1] = s.lo[2];
        base.shadow_box[2] = s.hi[0];
        base.shadow_box[3] = s.hi[2];
        base.shadow_plane[0] = s.lo[1] - 0.04f * (s.hi[1] - s.lo[1]);
        base.shadow_plane[1] = extent / 12.0f;
        base.shadow_plane[2] = 0.28f;
        if (!lighting.perspective) {
            std::array<double, 3> hv{};
            for (std::size_t i = 0; i < 3; ++i) hv[i] = double(base.light[i]) + double(lighting.eye[i]);
            hv = unitd(hv);
            for (std::size_t i = 0; i < 3; ++i) base.half_vector[i] = static_cast<float>(hv[i]);
        }
    }
    const Colormap map = s.colorbar ? s.colorbar->map : Colormap::turbo;
    for (int i = 0; i < 256; ++i) {
        const Rgb c = sample(map, i / 255.0);
        base.colormap[i][0] = c.r;
        base.colormap[i][1] = c.g;
        base.colormap[i][2] = c.b;
        base.colormap[i][3] = 1.0f;
    }
    base.no_data[0] = no_data.r;
    base.no_data[1] = no_data.g;
    base.no_data[2] = no_data.b;
    base.no_data[3] = 1.0f;
    std::vector<FrameConstants> constants;
    for (const MeshRange& r : ranges) {
        FrameConstants c = base;
        c.flags = r.flags;
        constants.push_back(c);
    }
    for (const LineRange& r : line_ranges) {
        FrameConstants c = base;
        c.line_width = r.width * supersample;
        constants.push_back(c);
    }
    if (constants.empty()) constants.push_back(base);

    auto cb = upload(constants.data(), constants.size() * sizeof(FrameConstants));
    if (!cb) return std::unexpected(std::move(cb.error()));

    const auto record_passes = [&]() {
        auto rtv = rtv_heap->GetCPUDescriptorHandleForHeapStart();
        const D3D12_CPU_DESCRIPTOR_HANDLE ssaa_rtv = rtv;
        const D3D12_CPU_DESCRIPTOR_HANDLE accum_rtv{rtv.ptr + rtv_size};
        const D3D12_CPU_DESCRIPTOR_HANDLE reveal_rtv{rtv.ptr + 2 * rtv_size};
        const D3D12_CPU_DESCRIPTOR_HANDLE final_rtv{rtv.ptr + 3 * rtv_size};
        const D3D12_CPU_DESCRIPTOR_HANDLE dsv = dsv_heap->GetCPUDescriptorHandleForHeapStart();
        const float white[4] = {1.0f, 1.0f, 1.0f, 1.0f};
        const float zero[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        const float one[4] = {1.0f, 1.0f, 1.0f, 1.0f};
        list->ClearRenderTargetView(ssaa_rtv, white, 0, nullptr);
        list->ClearRenderTargetView(accum_rtv, zero, 0, nullptr);
        list->ClearRenderTargetView(reveal_rtv, one, 0, nullptr);
        list->ClearDepthStencilView(dsv, D3D12_CLEAR_FLAG_DEPTH, 1.0f, 0, 0, nullptr);
        const D3D12_VIEWPORT vp{view.left * supersample, view.top * supersample, (view.right - view.left) * supersample,
                                (view.bottom - view.top) * supersample, 0.0f, 1.0f};
        const D3D12_RECT sc{static_cast<LONG>(vp.TopLeftX), static_cast<LONG>(vp.TopLeftY),
                            static_cast<LONG>(vp.TopLeftX + vp.Width), static_cast<LONG>(vp.TopLeftY + vp.Height)};
        list->SetGraphicsRootSignature(root.Get());
        ID3D12DescriptorHeap* heaps[] = {srv_heap.Get()};
        list->SetDescriptorHeaps(1, heaps);
        list->RSSetViewports(1, &vp);
        list->RSSetScissorRects(1, &sc);
        list->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        const D3D12_GPU_VIRTUAL_ADDRESS cb_base = (*cb)->GetGPUVirtualAddress();
        const D3D12_VERTEX_BUFFER_VIEW mesh_view{g.mesh->GetGPUVirtualAddress(),
                                                 static_cast<UINT>(std::max<std::size_t>(g.mesh_vertices * sizeof(MeshVertex), 256)),
                                                 sizeof(MeshVertex)};
        const D3D12_VERTEX_BUFFER_VIEW scene_line_view{g.lines->GetGPUVirtualAddress(),
                                                       static_cast<UINT>(std::max<std::size_t>(g.line_vertices * sizeof(LineVertex), 256)),
                                                       sizeof(LineVertex)};
        const D3D12_VERTEX_BUFFER_VIEW box_view{box_buffer ? box_buffer->GetGPUVirtualAddress() : 0,
                                                static_cast<UINT>(std::max<std::size_t>(box_lines.size() * sizeof(LineVertex), 256)),
                                                sizeof(LineVertex)};
        // A line range's buffer: the scene's, or the axis box's (Unit 26i).
        const auto line_view = [&](const LineRange& r) { return r.box ? &box_view : &scene_line_view; };
        if (enhanced) {
            list->OMSetRenderTargets(1, &ssaa_rtv, FALSE, nullptr);
            list->SetPipelineState(background.Get());
            list->SetGraphicsRootConstantBufferView(0, cb_base);
            list->DrawInstanced(3, 1, 0, 0);
            if (shadowed) {
                list->SetPipelineState(shadow.Get());
                list->DrawInstanced(6, 1, 0, 0);
            }
        }
        // Opaque triangles, then the depth-tested lines.
        list->OMSetRenderTargets(1, &ssaa_rtv, FALSE, &dsv);
        list->SetPipelineState(opaque.Get());
        list->IASetVertexBuffers(0, 1, &mesh_view);
        for (std::size_t k = 0; k < ranges.size(); ++k) {
            if (ranges[k].transparent || ranges[k].count == 0) continue;
            list->SetGraphicsRootConstantBufferView(0, cb_base + k * sizeof(FrameConstants));
            list->DrawInstanced(static_cast<UINT>(ranges[k].count), 1, static_cast<UINT>(ranges[k].first), 0);
        }
        list->SetPipelineState(line_depth.Get());
        for (std::size_t k = 0; k < line_ranges.size(); ++k) {
            if (line_ranges[k].top || line_ranges[k].count == 0) continue;
            list->IASetVertexBuffers(0, 1, line_view(line_ranges[k]));
            list->SetGraphicsRootConstantBufferView(0, cb_base + (ranges.size() + k) * sizeof(FrameConstants));
            list->DrawInstanced(static_cast<UINT>(line_ranges[k].count), 1, static_cast<UINT>(line_ranges[k].first), 0);
        }
        // Translucent triangles into the OIT targets, then composited over.
        bool any_transparent = false;
        for (const MeshRange& r : ranges) any_transparent = any_transparent || (r.transparent && r.count > 0);
        if (any_transparent) {
            const D3D12_CPU_DESCRIPTOR_HANDLE oit[] = {accum_rtv, reveal_rtv};
            list->OMSetRenderTargets(2, oit, FALSE, &dsv);
            list->SetPipelineState(transparent.Get());
            list->IASetVertexBuffers(0, 1, &mesh_view);
            for (std::size_t k = 0; k < ranges.size(); ++k) {
                if (!ranges[k].transparent || ranges[k].count == 0) continue;
                list->SetGraphicsRootConstantBufferView(0, cb_base + k * sizeof(FrameConstants));
                list->DrawInstanced(static_cast<UINT>(ranges[k].count), 1, static_cast<UINT>(ranges[k].first), 0);
            }
            D3D12_RESOURCE_BARRIER b[] = {
                transition(accum.Get(), D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE),
                transition(reveal.Get(), D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE)};
            list->ResourceBarrier(2, b);
            list->OMSetRenderTargets(1, &ssaa_rtv, FALSE, nullptr);
            list->SetPipelineState(composite.Get());
            list->SetGraphicsRootConstantBufferView(0, cb_base);
            list->SetGraphicsRootDescriptorTable(1, srv_heap->GetGPUDescriptorHandleForHeapStart());
            const D3D12_VIEWPORT full{0.0f, 0.0f, static_cast<float>(width * supersample), static_cast<float>(height * supersample), 0.0f, 1.0f};
            const D3D12_RECT full_sc{0, 0, width * supersample, height * supersample};
            list->RSSetViewports(1, &full);
            list->RSSetScissorRects(1, &full_sc);
            list->DrawInstanced(3, 1, 0, 0);
            D3D12_RESOURCE_BARRIER a[] = {
                transition(accum.Get(), D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE, D3D12_RESOURCE_STATE_RENDER_TARGET),
                transition(reveal.Get(), D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE, D3D12_RESOURCE_STATE_RENDER_TARGET)};
            list->ResourceBarrier(2, a);
            list->RSSetViewports(1, &vp);
            list->RSSetScissorRects(1, &sc);
        }
        // Lines drawn over everything (2D outlines, contours).
        list->OMSetRenderTargets(1, &ssaa_rtv, FALSE, &dsv);
        list->SetPipelineState(line_top.Get());
        for (std::size_t k = 0; k < line_ranges.size(); ++k) {
            if (!line_ranges[k].top || line_ranges[k].count == 0) continue;
            list->IASetVertexBuffers(0, 1, line_view(line_ranges[k]));
            list->SetGraphicsRootConstantBufferView(0, cb_base + (ranges.size() + k) * sizeof(FrameConstants));
            list->DrawInstanced(static_cast<UINT>(line_ranges[k].count), 1, static_cast<UINT>(line_ranges[k].first), 0);
        }
        // Average 2 x 2 into the output.
        {
            const auto b = transition(ssaa.Get(), D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
            list->ResourceBarrier(1, &b);
            if (final_state != D3D12_RESOURCE_STATE_RENDER_TARGET) {
                const auto f = transition(final_target.Get(), final_state, D3D12_RESOURCE_STATE_RENDER_TARGET);
                list->ResourceBarrier(1, &f);
                final_state = D3D12_RESOURCE_STATE_RENDER_TARGET;
            }
            list->OMSetRenderTargets(1, &final_rtv, FALSE, nullptr);
            list->SetPipelineState(resolve.Get());
            list->SetGraphicsRootConstantBufferView(0, cb_base);
            D3D12_GPU_DESCRIPTOR_HANDLE table = srv_heap->GetGPUDescriptorHandleForHeapStart();
            table.ptr += 2 * srv_size;
            list->SetGraphicsRootDescriptorTable(1, table);
            const D3D12_VIEWPORT out{0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height), 0.0f, 1.0f};
            const D3D12_RECT out_sc{0, 0, width, height};
            list->RSSetViewports(1, &out);
            list->RSSetScissorRects(1, &out_sc);
            list->DrawInstanced(3, 1, 0, 0);
            const auto back = transition(ssaa.Get(), D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE, D3D12_RESOURCE_STATE_RENDER_TARGET);
            list->ResourceBarrier(1, &back);
        }
    };
    record_passes();
    if (auto r = execute(); !r) return std::unexpected(std::move(r.error()));
    // frame_time: the passes again on the resident geometry, what redrawing a window costs.
    if (repeat_frames > 0) {
        const auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < repeat_frames; ++k) {
            record_passes();
            if (auto r = execute(); !r) return std::unexpected(std::move(r.error()));
        }
        repeat_ms = 1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / repeat_frames;
    }

    // Overlays on Direct2D.
    if (auto r = begin_2d(); !r) return std::unexpected(std::move(r.error()));
    {
        Ink ink(d2d.Get(), dwrite.Get(), &layouts);
        d2d->SetAntialiasMode(D2D1_ANTIALIAS_MODE_PER_PRIMITIVE);
        if (!s.title.empty()) ink.text(s.title, title_text.Get(), 14.0f, 10.0f, 0, 0);
        if (s.colorbar) {
            draw_colorbar(ink, text.Get(), small_text.Get(), *s.colorbar,
                          {static_cast<float>(w) - bar_w + 8.0f, title_h, static_cast<float>(w) - 8.0f,
                           static_cast<float>(h) - 12.0f}, enhanced);
        }
        const D2D1_COLOR_F axis{0.25f, 0.27f, 0.3f, 1.0f};
        if (axes) {
            // The x and y axes of the 2D map in um, along the view's bottom and left.
            const D2D1_POINT_2F p0 = project(v, s.lo, view), p1 = project(v, s.hi, view);
            float xl = std::min(p0.x, p1.x), xr = std::max(p0.x, p1.x);
            float yb = std::max(p0.y, p1.y), yt = std::min(p0.y, p1.y);
            // Zoomed or panned (Unit 26): the axes span the part of the bounds in view. A fitted view
            // holds the whole bounds, and these are the bounds' own values.
            double x_lo = s.lo[0], x_hi = s.hi[0], y_lo = s.lo[1], y_hi = s.hi[1];
            if (xl < view.left || xr > view.right) {
                const auto world_x = [&](float px) { return s.lo[0] + (px - p0.x) * (s.hi[0] - s.lo[0]) / (p1.x - p0.x); };
                x_lo = std::max(x_lo, static_cast<double>(world_x(view.left)));
                x_hi = std::min(x_hi, static_cast<double>(world_x(view.right)));
                xl = std::max(xl, view.left);
                xr = std::min(xr, view.right);
            }
            if (yt < view.top || yb > view.bottom) {
                const auto world_y = [&](float py) { return s.lo[1] + (py - p0.y) * (s.hi[1] - s.lo[1]) / (p1.y - p0.y); };
                y_lo = std::max(y_lo, static_cast<double>(world_y(view.bottom)));
                y_hi = std::min(y_hi, static_cast<double>(world_y(view.top)));
                yb = std::min(yb, view.bottom);
                yt = std::max(yt, view.top);
            }
            if (xl < xr) d2d->DrawLine({xl, yb + 6.0f}, {xr, yb + 6.0f}, ink.color(axis), 1.0f);
            if (yt < yb) d2d->DrawLine({xl - 6.0f, yt}, {xl - 6.0f, yb}, ink.color(axis), 1.0f);
            for (const Tick& t : x_lo < x_hi ? linear_ticks(x_lo, x_hi, 6) : std::vector<Tick>{}) {
                const float x = project(v, {static_cast<float>(t.value), s.lo[1], 0.0f}, view).x;
                d2d->DrawLine({x, yb + 6.0f}, {x, yb + 11.0f}, ink.color(axis), 1.0f);
                ink.text(t.label, small_text.Get(), x, yb + 13.0f, 1, 0);
            }
            for (const Tick& t : y_lo < y_hi ? linear_ticks(y_lo, y_hi, 5) : std::vector<Tick>{}) {
                const float y = project(v, {s.lo[0], static_cast<float>(t.value), 0.0f}, view).y;
                d2d->DrawLine({xl - 11.0f, y}, {xl - 6.0f, y}, ink.color(axis), 1.0f);
                ink.text(t.label, small_text.Get(), xl - 14.0f, y, 2, 1);
            }
            ink.text("x (\xC2\xB5m)", text.Get(), 0.5f * (xl + xr), yb + 30.0f, 1, 0);
            d2d->SetTransform(D2D1::Matrix3x2F::Rotation(-90.0f, {14.0f, 0.5f * (yt + yb)}));
            ink.text("y (\xC2\xB5m)", text.Get(), 14.0f, 0.5f * (yt + yb), 1, 1);
            d2d->SetTransform(D2D1::Matrix3x2F::Identity());
        }
        // Unit 26h (JD1): the axis box's ticks, labels and titles, outside the device, inside the view.
        const auto inside = [&](float x, float y) {
            return x >= view.left && x <= view.right && y >= view.top && y <= view.bottom;
        };
        const D2D1_COLOR_F box_ink{0.35f, 0.37f, 0.41f, 1.0f};
        for (const BoxAxis& ax : box.axes) {
            for (const AxisLabel& l : ax.ticks) {
                if (!inside(l.x, l.y)) continue;
                d2d->DrawLine({l.x, l.y}, {l.x + ax.outward[0] * 6.0f, l.y + ax.outward[1] * 6.0f}, ink.color(box_ink), 1.0f);
                if (l.kept && inside(l.label_x, l.label_y)) ink.text(l.text, small_text.Get(), l.label_x, l.label_y, 1, 1, box_ink);
            }
            if (inside(ax.title_x, ax.title_y)) ink.text(ax.title, small_text.Get(), ax.title_x, ax.title_y, 1, 1, box_ink);
        }
        // Unit 26h (JD2): the view cube in place of the triad in the enhanced look.
        const ViewCube cube = enhanced ? view_cube(s, w, h) : ViewCube{};
        if (!cube.faces.empty()) {
            for (const CubeFace& f : cube.faces) {
                ComPtr<ID2D1PathGeometry> path;
                ComPtr<ID2D1GeometrySink> sink;
                if (FAILED(d2d_factory->CreatePathGeometry(&path)) || FAILED(path->Open(&sink))) continue;
                sink->BeginFigure({f.polygon[0][0], f.polygon[0][1]}, D2D1_FIGURE_BEGIN_FILLED);
                for (std::size_t k = 1; k < 4; ++k) sink->AddLine({f.polygon[k][0], f.polygon[k][1]});
                sink->EndFigure(D2D1_FIGURE_END_CLOSED);
                if (FAILED(sink->Close())) continue;
                const D2D1_COLOR_F fill{0.90f * f.shade, 0.92f * f.shade, 0.95f * f.shade, 1.0f};
                d2d->FillGeometry(path.Get(), ink.color(fill));
                d2d->DrawGeometry(path.Get(), ink.color({0.42f, 0.45f, 0.5f, 1.0f}), 1.0f);
                // The name where the face is wide and tall enough to hold it (6 px a character).
                float x0 = 1e9f, x1 = -1e9f, y0 = 1e9f, y1 = -1e9f;
                for (const auto& p : f.polygon) {
                    x0 = std::min(x0, p[0]);
                    x1 = std::max(x1, p[0]);
                    y0 = std::min(y0, p[1]);
                    y1 = std::max(y1, p[1]);
                }
                if (x1 - x0 >= 6.0f * static_cast<float>(f.label.size()) + 4.0f && y1 - y0 >= 16.0f) {
                    float cx = 0.0f, cy = 0.0f;
                    for (const auto& p : f.polygon) {
                        cx += 0.25f * p[0];
                        cy += 0.25f * p[1];
                    }
                    ink.text(f.label, small_text.Get(), cx, cy, 1, 1, {0.2f, 0.22f, 0.26f, 1.0f});
                }
            }
            const std::array<D2D1_COLOR_F, 3> colors{
                {{0.80f, 0.20f, 0.18f, 1.0f}, {0.18f, 0.60f, 0.25f, 1.0f}, {0.20f, 0.36f, 0.80f, 1.0f}}};
            for (std::size_t a = 0; a < cube.arrows.size(); ++a) {
                const CubeArrow& r = cube.arrows[a];
                if (r.from == r.to) continue;
                d2d->DrawLine({r.from[0], r.from[1]}, {r.to[0], r.to[1]}, ink.color(colors[a]), 2.0f);
                const float dx = r.to[0] - r.from[0], dy = r.to[1] - r.from[1], n = std::hypot(dx, dy);
                ink.text(r.name, small_text.Get(), r.to[0] + 8.0f * dx / n, r.to[1] + 8.0f * dy / n, 1, 1, colors[a]);
            }
        } else if (!axes && s.axis_triad) {
            const D2D1_POINT_2F o{view.left + 44.0f, view.bottom - 44.0f};
            const std::array<std::pair<const char*, D2D1_COLOR_F>, 3> axes3{
                {{"x", {0.80f, 0.20f, 0.18f, 1.0f}}, {"y", {0.18f, 0.60f, 0.25f, 1.0f}}, {"z", {0.20f, 0.36f, 0.80f, 1.0f}}}};
            for (int a = 0; a < 3; ++a) {
                std::array<double, 3> e{};
                e[static_cast<std::size_t>(a)] = 1.0;
                const float dx = static_cast<float>(dotd(e, v.right)) * 30.0f;
                const float dy = -static_cast<float>(dotd(e, v.up)) * 30.0f;
                const auto& [name, c] = axes3[static_cast<std::size_t>(a)];
                d2d->DrawLine(o, {o.x + dx, o.y + dy}, ink.color(c), 2.0f);
                ink.text(name, small_text.Get(), o.x + dx * 1.3f, o.y + dy * 1.3f, 1, 1, c);
            }
        }
        // Unit 26o (BSD6): the band surface's plane in its device, in the top right corner.
        if (const InsetFigure inset = plane_inset_figure(s, w, h); !inset.edges.empty()) {
            const auto& r = inset.frame;
            d2d->FillRoundedRectangle({{r[0], r[1], r[2], r[3]}, 6.0f, 6.0f}, ink.color({1.0f, 1.0f, 1.0f, 0.82f}));
            d2d->DrawRoundedRectangle({{r[0], r[1], r[2], r[3]}, 6.0f, 6.0f}, ink.color({0.62f, 0.65f, 0.7f, 1.0f}), 1.0f);
            ComPtr<ID2D1PathGeometry> path;
            ComPtr<ID2D1GeometrySink> sink;
            if (SUCCEEDED(d2d_factory->CreatePathGeometry(&path)) && SUCCEEDED(path->Open(&sink))) {
                sink->BeginFigure({inset.plane[0][0], inset.plane[0][1]}, D2D1_FIGURE_BEGIN_FILLED);
                for (std::size_t k = 1; k < 4; ++k) sink->AddLine({inset.plane[k][0], inset.plane[k][1]});
                sink->EndFigure(D2D1_FIGURE_END_CLOSED);
                if (SUCCEEDED(sink->Close())) {
                    d2d->FillGeometry(path.Get(), ink.color({0.95f, 0.55f, 0.15f, 0.40f}));
                    d2d->DrawGeometry(path.Get(), ink.color({0.80f, 0.38f, 0.05f, 1.0f}), 1.5f);
                }
            }
            for (const auto& e : inset.edges) {
                d2d->DrawLine({e[0][0], e[0][1]}, {e[1][0], e[1][1]}, ink.color({0.30f, 0.33f, 0.38f, 1.0f}), 1.0f);
            }
        }
        // Unit 26n (DT2): the scene's labels, centred just above their points, over everything.
        for (const SceneLabel& l : s.labels) {
            const D2D1_POINT_2F q = project(v, l.at, view);
            ink.text(l.text, text.Get(), q.x, q.y - 3.0f, 1, 2, {l.color.r, l.color.g, l.color.b, l.color.a});
        }
    }
    return end_2d();
}

PlotFrame plot_frame(const PlotScene& s, int w, int h) noexcept {
    const float top = s.title.empty() ? 18.0f : 44.0f;
    return {86.0f, top, static_cast<float>(w) - 24.0f, static_cast<float>(h) - 58.0f};
}

namespace {

// A value's place along an axis, 0 at lo and 1 at hi (NaN off a log axis).
double axis_fraction(const Axis& a, double v) noexcept {
    if (a.scale == AxisScale::log10) {
        if (!(v > 0.0)) return std::numeric_limits<double>::quiet_NaN();
        return (std::log10(v) - std::log10(a.lo)) / (std::log10(a.hi) - std::log10(a.lo));
    }
    return (v - a.lo) / (a.hi - a.lo);
}
float plot_x(const PlotFrame& f, const Axis& a, double x) noexcept {
    return f.left + static_cast<float>(axis_fraction(a, x)) * (f.right - f.left);
}
float plot_y(const PlotFrame& f, const Axis& a, double y) noexcept {
    return f.bottom - static_cast<float>(axis_fraction(a, y)) * (f.bottom - f.top);
}

}  // namespace

ScreenPoint plot_position(const PlotScene& s, int w, int h, double x, double y) noexcept {
    const PlotFrame f = plot_frame(s, w, h);
    const float sx = plot_x(f, s.x, x), sy = plot_y(f, s.y, y);
    if (!std::isfinite(sx) || !std::isfinite(sy)) {
        const float nan = std::numeric_limits<float>::quiet_NaN();
        return {nan, nan, 0.0f};
    }
    return {sx, sy, 0.0f};
}

std::expected<void, base::Error> Renderer::Impl::render_plot(const PlotScene& s, int w, int h) {
    if (auto r = ensure_targets(w, h); !r) return std::unexpected(std::move(r.error()));
    if (final_state != D3D12_RESOURCE_STATE_RENDER_TARGET) {
        const auto f = transition(final_target.Get(), final_state, D3D12_RESOURCE_STATE_RENDER_TARGET);
        list->ResourceBarrier(1, &f);
        final_state = D3D12_RESOURCE_STATE_RENDER_TARGET;
    }
    const float white[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    auto rtv = rtv_heap->GetCPUDescriptorHandleForHeapStart();
    rtv.ptr += 3 * rtv_size;
    list->ClearRenderTargetView(rtv, white, 0, nullptr);
    if (auto r = execute(); !r) return std::unexpected(std::move(r.error()));
    if (auto r = begin_2d(); !r) return std::unexpected(std::move(r.error()));
    {
        Ink ink(d2d.Get(), dwrite.Get(), &layouts);
        d2d->SetAntialiasMode(D2D1_ANTIALIAS_MODE_PER_PRIMITIVE);
        const PlotFrame f = plot_frame(s, w, h);
        const D2D1_RECT_F frame{f.left, f.top, f.right, f.bottom};
        const auto px = [&](double x) { return plot_x(f, s.x, x); };
        const auto py = [&](double y) { return plot_y(f, s.y, y); };
        const D2D1_COLOR_F grid{0.89f, 0.91f, 0.93f, 1.0f}, axis{0.25f, 0.27f, 0.3f, 1.0f};
        for (const Tick& t : s.x.ticks) {
            const float x = frame.left + static_cast<float>(t.position) * (frame.right - frame.left);
            if (t.major) d2d->DrawLine({x, frame.top}, {x, frame.bottom}, ink.color(grid), 1.0f);
            d2d->DrawLine({x, frame.bottom}, {x, frame.bottom + (t.major ? 5.0f : 3.0f)}, ink.color(axis), 1.0f);
            if (t.major) ink.text(t.label, small_text.Get(), x, frame.bottom + 8.0f, 1, 0);
        }
        for (const Tick& t : s.y.ticks) {
            const float y = frame.bottom - static_cast<float>(t.position) * (frame.bottom - frame.top);
            if (t.major) d2d->DrawLine({frame.left, y}, {frame.right, y}, ink.color(grid), 1.0f);
            d2d->DrawLine({frame.left - (t.major ? 5.0f : 3.0f), y}, {frame.left, y}, ink.color(axis), 1.0f);
            if (t.major) ink.text(t.label, small_text.Get(), frame.left - 8.0f, y, 2, 1);
        }
        d2d->DrawRectangle(frame, ink.color(axis), 1.0f);
        d2d->PushAxisAlignedClip(frame, D2D1_ANTIALIAS_MODE_ALIASED);
        for (const Series& series : s.series) {
            const D2D1_COLOR_F c = d2(series.color);
            for (std::size_t i = 0; i + 1 < series.x.size(); ++i) {
                const float x0 = px(series.x[i]), y0 = py(series.y[i]);
                const float x1 = px(series.x[i + 1]), y1 = py(series.y[i + 1]);
                if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(x1) || !std::isfinite(y1)) continue;
                d2d->DrawLine({x0, y0}, {x1, y1}, ink.color(c), 2.0f);
            }
            if (series.markers) {
                for (std::size_t i = 0; i < series.x.size(); ++i) {
                    const float x = px(series.x[i]), y = py(series.y[i]);
                    if (!std::isfinite(x) || !std::isfinite(y)) continue;
                    const D2D1_ELLIPSE e{{x, y}, 3.5f, 3.5f};
                    const bool open = i < series.flagged.size() && series.flagged[i];
                    d2d->FillEllipse(e, ink.color(open ? D2D1_COLOR_F{1, 1, 1, 1} : c));
                    if (open) d2d->DrawEllipse(e, ink.color(c), 1.5f);
                }
            }
        }
        d2d->PopAxisAlignedClip();
        ink.text(s.x.label, text.Get(), 0.5f * (frame.left + frame.right), static_cast<float>(h) - 8.0f, 1, 2);
        d2d->SetTransform(D2D1::Matrix3x2F::Rotation(-90.0f, {16.0f, 0.5f * (frame.top + frame.bottom)}));
        ink.text(s.y.label, text.Get(), 16.0f, 0.5f * (frame.top + frame.bottom), 1, 1);
        d2d->SetTransform(D2D1::Matrix3x2F::Identity());
        if (!s.title.empty()) ink.text(s.title, title_text.Get(), frame.left, 12.0f, 0, 0);
        if (s.legend && !s.series.empty()) {
            float wmax = 0.0f;
            for (const Series& series : s.series) wmax = std::max(wmax, ink.width_of(series.name, small_text.Get()));
            const float lw = wmax + 44.0f, lh = 8.0f + 18.0f * static_cast<float>(s.series.size());
            const D2D1_RECT_F box{frame.right - lw - 8.0f, frame.top + 8.0f, frame.right - 8.0f, frame.top + 8.0f + lh};
            d2d->FillRectangle(box, ink.color({1, 1, 1, 0.92f}));
            d2d->DrawRectangle(box, ink.color(grid), 1.0f);
            for (std::size_t k = 0; k < s.series.size(); ++k) {
                const float y = box.top + 13.0f + 18.0f * static_cast<float>(k);
                d2d->DrawLine({box.left + 8.0f, y}, {box.left + 30.0f, y}, ink.color(d2(s.series[k].color)), 2.0f);
                ink.text(s.series[k].name, small_text.Get(), box.left + 36.0f, y, 0, 1);
            }
        }
    }
    return end_2d();
}

// --- Renderer -----------------------------------------------------------------------------------

Renderer::Renderer(std::unique_ptr<Impl> impl) noexcept : impl_(std::move(impl)) {}
Renderer::Renderer(Renderer&&) noexcept = default;
Renderer& Renderer::operator=(Renderer&&) noexcept = default;
Renderer::~Renderer() = default;

std::expected<Renderer, base::Error> Renderer::create(const RendererOptions& options) {
    try {
        auto impl = std::make_unique<Impl>();
        if (auto r = impl->create(options); !r) return std::unexpected(std::move(r.error()));
        return Renderer(std::move(impl));
    } catch (const std::bad_alloc&) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "out of memory", std::nullopt});
    }
}

namespace {

std::expected<void, base::Error> check_size(int w, int h) {
    if (w < 16 || h < 16 || w > 16384 || h > 16384) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, "an image must be 16 to 16384 pixels each way",
                                           base::ErrorContext{.index = std::nullopt, .value = static_cast<double>(std::min(w, h))}});
    }
    return {};
}

}  // namespace

std::expected<Image, base::Error> Renderer::draw(const FieldScene& scene, int width, int height) {
    if (auto r = check_size(width, height); !r) return std::unexpected(std::move(r.error()));
    try {
        if (auto r = impl_->render_field(scene, width, height); !r) return std::unexpected(std::move(r.error()));
        return impl_->read_back();
    } catch (const std::bad_alloc&) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "out of memory", std::nullopt});
    }
}

std::expected<Image, base::Error> Renderer::draw(const PlotScene& scene, int width, int height) {
    if (auto r = check_size(width, height); !r) return std::unexpected(std::move(r.error()));
    try {
        if (auto r = impl_->render_plot(scene, width, height); !r) return std::unexpected(std::move(r.error()));
        return impl_->read_back();
    } catch (const std::bad_alloc&) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "out of memory", std::nullopt});
    }
}

std::expected<double, base::Error> Renderer::frame_time(const FieldScene& scene, int width, int height,
                                                        int frames) {
    if (auto r = check_size(width, height); !r) return std::unexpected(std::move(r.error()));
    impl_->repeat_frames = std::max(frames, 1);
    impl_->repeat_ms = 0.0;
    auto drawn = draw(scene, width, height);
    impl_->repeat_frames = 0;
    if (!drawn) return std::unexpected(std::move(drawn.error()));
    return impl_->repeat_ms;
}

const std::string& Renderer::adapter_name() const noexcept { return impl_->adapter_name; }
bool Renderer::is_warp() const noexcept { return impl_->warp; }

std::size_t Renderer::debug_errors() const { return debug_messages().size(); }

std::vector<std::string> Renderer::debug_messages() const {
    std::vector<std::string> out;
    if (!impl_->info) return out;
    const UINT64 n = impl_->info->GetNumStoredMessages();
    for (UINT64 i = 0; i < n; ++i) {
        SIZE_T size = 0;
        impl_->info->GetMessage(i, nullptr, &size);
        std::vector<std::byte> buffer(size);
        auto* m = reinterpret_cast<D3D12_MESSAGE*>(buffer.data());
        if (SUCCEEDED(impl_->info->GetMessage(i, m, &size)) &&
            (m->Severity == D3D12_MESSAGE_SEVERITY_ERROR || m->Severity == D3D12_MESSAGE_SEVERITY_CORRUPTION)) {
            out.push_back(std::format("[{}] ", static_cast<int>(m->ID)) +
                          std::string(m->pDescription, m->DescriptionByteLength > 0 ? m->DescriptionByteLength - 1 : 0));
        }
    }
    return out;
}

// --- windows (Unit 26) ---------------------------------------------------------------------------

WindowTarget::WindowTarget(std::unique_ptr<Impl> impl) noexcept : impl_(std::move(impl)) {}
WindowTarget::WindowTarget(WindowTarget&&) noexcept = default;
WindowTarget& WindowTarget::operator=(WindowTarget&&) noexcept = default;
WindowTarget::~WindowTarget() = default;
int WindowTarget::width() const noexcept { return impl_->width; }
int WindowTarget::height() const noexcept { return impl_->height; }

std::expected<void, base::Error> Renderer::Impl::show(WindowTarget::Impl& t) {
    // The drawn image into the back buffer, then presented: the window shows what draw() reads back.
    ID3D12Resource* back = t.buffers[t.swap->GetCurrentBackBufferIndex()].Get();
    if (final_state != D3D12_RESOURCE_STATE_COPY_SOURCE) {
        const auto b = transition(final_target.Get(), final_state, D3D12_RESOURCE_STATE_COPY_SOURCE);
        list->ResourceBarrier(1, &b);
    }
    const auto in = transition(back, D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_COPY_DEST);
    list->ResourceBarrier(1, &in);
    list->CopyResource(back, final_target.Get());
    const D3D12_RESOURCE_BARRIER out[] = {
        transition(back, D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_PRESENT),
        transition(final_target.Get(), D3D12_RESOURCE_STATE_COPY_SOURCE, D3D12_RESOURCE_STATE_RENDER_TARGET)};
    list->ResourceBarrier(2, out);
    final_state = D3D12_RESOURCE_STATE_RENDER_TARGET;
    if (auto r = execute(); !r) return r;
    const HRESULT hr = t.swap->Present(1, 0);
    if (FAILED(hr)) return std::unexpected(device_error("IDXGISwapChain::Present", hr));
    return {};
}

std::expected<WindowTarget, base::Error> Renderer::window_target(void* window, int width, int height) {
    if (auto r = check_size(width, height); !r) return std::unexpected(std::move(r.error()));
    auto t = std::make_unique<WindowTarget::Impl>();
    t->window = static_cast<HWND>(window);
    t->flush = [impl = impl_.get()] { (void)impl->wait(); };
    DXGI_SWAP_CHAIN_DESC1 d{};
    d.Width = static_cast<UINT>(width);
    d.Height = static_cast<UINT>(height);
    d.Format = color_format;
    d.SampleDesc.Count = 1;
    d.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    d.BufferCount = 2;
    d.Scaling = DXGI_SCALING_NONE;
    d.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
    d.AlphaMode = DXGI_ALPHA_MODE_IGNORE;
    ComPtr<IDXGISwapChain1> s1;
    NITCAD_HR(impl_->factory->CreateSwapChainForHwnd(impl_->queue.Get(), t->window, &d, nullptr, nullptr, &s1),
              "IDXGIFactory2::CreateSwapChainForHwnd");
    NITCAD_HR(s1.As(&t->swap), "QueryInterface IDXGISwapChain3");
    impl_->factory->MakeWindowAssociation(t->window, DXGI_MWA_NO_ALT_ENTER);
    for (UINT i = 0; i < 2; ++i) NITCAD_HR(t->swap->GetBuffer(i, IID_PPV_ARGS(&t->buffers[i])), "IDXGISwapChain::GetBuffer");
    t->width = width;
    t->height = height;
    return WindowTarget(std::move(t));
}

namespace {

// The swap chain at the size: its buffers released, resized and taken again.
std::expected<void, base::Error> fit_target(WindowTarget::Impl& t, int width, int height) {
    if (t.width == width && t.height == height) return {};
    for (auto& b : t.buffers) b.Reset();
    NITCAD_HR(t.swap->ResizeBuffers(2, static_cast<UINT>(width), static_cast<UINT>(height), color_format, 0),
              "IDXGISwapChain::ResizeBuffers");
    for (UINT i = 0; i < 2; ++i) NITCAD_HR(t.swap->GetBuffer(i, IID_PPV_ARGS(&t.buffers[i])), "IDXGISwapChain::GetBuffer");
    t.width = width;
    t.height = height;
    return {};
}

}  // namespace

std::expected<void, base::Error> Renderer::present(const FieldScene& scene, WindowTarget& target, int width, int height) {
    if (auto r = check_size(width, height); !r) return std::unexpected(std::move(r.error()));
    try {
        if (auto r = impl_->wait(); !r) return r;  // the swap chain's buffers are free
        if (auto r = fit_target(*target.impl_, width, height); !r) return r;
        if (auto r = impl_->render_field(scene, width, height, &target.impl_->resident); !r) return r;
        return impl_->show(*target.impl_);
    } catch (const std::bad_alloc&) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "out of memory", std::nullopt});
    }
}

std::expected<void, base::Error> Renderer::present_camera(const FieldScene& scene, WindowTarget& target, int width,
                                                          int height) {
    if (auto r = check_size(width, height); !r) return std::unexpected(std::move(r.error()));
    try {
        if (auto r = impl_->wait(); !r) return r;
        if (auto r = fit_target(*target.impl_, width, height); !r) return r;
        if (auto r = impl_->render_field(scene, width, height, &target.impl_->resident, true); !r) return r;
        return impl_->show(*target.impl_);
    } catch (const std::bad_alloc&) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "out of memory", std::nullopt});
    }
}

bool Renderer::last_present_reused() const noexcept { return impl_->reused; }

std::expected<Image, base::Error> Renderer::last_image() {
    if (impl_->width == 0) return std::unexpected(base::Error{base::ErrorCode::invalid_input, "nothing drawn yet", std::nullopt});
    return impl_->read_back();
}

std::expected<void, base::Error> Renderer::present(const PlotScene& scene, WindowTarget& target, int width, int height) {
    if (auto r = check_size(width, height); !r) return std::unexpected(std::move(r.error()));
    try {
        if (auto r = impl_->wait(); !r) return r;
        if (auto r = fit_target(*target.impl_, width, height); !r) return r;
        target.impl_->resident = {};  // a plot: the field scene's buffers are no longer drawn (NG5)
        if (auto r = impl_->render_plot(scene, width, height); !r) return r;
        return impl_->show(*target.impl_);
    } catch (const std::bad_alloc&) {
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted, "out of memory", std::nullopt});
    }
}

ScreenPoint screen_position(const FieldScene& scene, int width, int height, const Float3& p) {
    const D2D1_RECT_F view = field_view_rect(scene, width, height);
    const View v = make_view(scene, (view.right - view.left) / (view.bottom - view.top));
    const Mat4& m = v.view_projection;
    const float z = p[0] * m[2] + p[1] * m[6] + p[2] * m[10] + m[14];
    const float w = p[0] * m[3] + p[1] * m[7] + p[2] * m[11] + m[15];
    const D2D1_POINT_2F px = project(v, p, view);
    return {px.x, px.y, w != 0.0f ? z / w : z};
}

// --- the view in double precision (Unit 26l, SD4) ------------------------------------------------------------------

// make_view's layout, kept in double: the same expressions, without the float matrix.
ViewFrame view_frame(const FieldScene& s, int width, int height) {
    const D2D1_RECT_F rect = field_view_rect(s, width, height);
    const double aspect = (rect.right - rect.left) / (rect.bottom - rect.top);
    const double az = s.camera.azimuth_deg * std::numbers::pi / 180.0;
    const double el = s.camera.elevation_deg * std::numbers::pi / 180.0;
    const std::array<double, 3> eye{std::cos(el) * std::sin(az), std::sin(el), std::cos(el) * std::cos(az)};
    ViewFrame f;
    f.forward = {-eye[0], -eye[1], -eye[2]};
    std::array<double, 3> world_up{0.0, 1.0, 0.0};
    if (std::abs(dotd(world_up, f.forward)) > 0.999) world_up = {0.0, 0.0, -1.0};
    f.right = unitd(crossd(f.forward, world_up));
    f.up = crossd(f.right, f.forward);
    const std::array<double, 3> lo{s.lo[0], s.lo[1], s.lo[2]}, hi{s.hi[0], s.hi[1], s.hi[2]};
    f.centre = {0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300, radius = 0.0;
    for (int k = 0; k < 8; ++k) {
        const std::array<double, 3> p{(k & 1) ? hi[0] : lo[0], (k & 2) ? hi[1] : lo[1], (k & 4) ? hi[2] : lo[2]};
        const std::array<double, 3> d{p[0] - f.centre[0], p[1] - f.centre[1], p[2] - f.centre[2]};
        x0 = std::min(x0, dotd(d, f.right));
        x1 = std::max(x1, dotd(d, f.right));
        y0 = std::min(y0, dotd(d, f.up));
        y1 = std::max(y1, dotd(d, f.up));
        radius = std::max(radius, std::sqrt(dotd(d, d)));
    }
    f.radius = std::max(radius, 1e-6);
    f.half_h = std::max({0.5 * (y1 - y0), 0.5 * (x1 - x0) / aspect, 1e-6}) / 0.92 / static_cast<double>(s.camera.zoom);
    f.half_w = f.half_h * aspect;
    f.cx = 0.5 * (x0 + x1) - static_cast<double>(s.camera.pan[0]) * f.half_w;
    f.cy = 0.5 * (y0 + y1) - static_cast<double>(s.camera.pan[1]) * f.half_h;
    f.left = rect.left;
    f.top = rect.top;
    f.right_px = rect.right;
    f.bottom_px = rect.bottom;
    if (s.camera.projection == Projection::perspective) {
        f.perspective = true;
        f.tan_half = std::tan(0.5 * s.camera.fov_deg * std::numbers::pi / 180.0);
        f.depth = f.half_h / f.tan_half + f.radius;
        for (std::size_t i = 0; i < 3; ++i) f.eye[i] = f.centre[i] + f.cx * f.right[i] + f.cy * f.up[i] - f.depth * f.forward[i];
    }
    return f;
}

ScreenRay screen_ray(const FieldScene& s, int width, int height, double x, double y) {
    const ViewFrame f = view_frame(s, width, height);
    const double nx = 2.0 * (x - f.left) / (f.right_px - f.left) - 1.0;
    const double ny = 1.0 - 2.0 * (y - f.top) / (f.bottom_px - f.top);
    ScreenRay r;
    if (!f.perspective) {
        // vx = nx half_w, vy = ny half_h; p = c + (vx + cx) right + (vy + cy) up + vz forward, from in front.
        const double vx = nx * f.half_w + f.cx, vy = ny * f.half_h + f.cy, vz = -2.0 * f.radius;
        for (std::size_t i = 0; i < 3; ++i) r.origin[i] = f.centre[i] + vx * f.right[i] + vy * f.up[i] + vz * f.forward[i];
        r.direction = f.forward;
        return r;
    }
    // On the view plane through the centre (vz = 0): vx = nx D tan aspect, vy = ny D tan.
    const double aspect = f.half_w / f.half_h;
    const double vx = nx * f.depth * f.tan_half * aspect + f.cx, vy = ny * f.depth * f.tan_half + f.cy;
    std::array<double, 3> p{};
    for (std::size_t i = 0; i < 3; ++i) p[i] = f.centre[i] + vx * f.right[i] + vy * f.up[i];
    r.origin = f.eye;
    r.direction = unitd({p[0] - f.eye[0], p[1] - f.eye[1], p[2] - f.eye[2]});
    return r;
}

double um_per_pixel(const FieldScene& s, int width, int height, const std::array<double, 3>& p) {
    const ViewFrame f = view_frame(s, width, height);
    const double per_px = 2.0 * f.half_h / (f.bottom_px - f.top);  // on the view plane through the centre
    if (!f.perspective) return per_px;
    const double d = dotd({p[0] - f.eye[0], p[1] - f.eye[1], p[2] - f.eye[2]}, f.forward);
    return per_px * std::max(d, 1e-9) / f.depth;
}

// --- the axis box (Unit 26h, JD1) -----------------------------------------------------------------

AxisBox axis_box(const FieldScene& s, int width, int height) {
    AxisBox out;
    if (s.dimension != 3) return out;
    const D2D1_RECT_F view = field_view_rect(s, width, height);
    const View v = make_view(s, (view.right - view.left) / (view.bottom - view.top));
    const std::array<double, 3> lo{s.lo[0], s.lo[1], s.lo[2]}, hi{s.hi[0], s.hi[1], s.hi[2]};
    const double extent = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-30});
    // Toward the eye from a point: the fixed direction (orthographic) or the camera's position.
    const auto toward_eye = [&](const std::array<double, 3>& p) {
        if (!v.perspective) return std::array<double, 3>{-v.forward[0], -v.forward[1], -v.forward[2]};
        return std::array<double, 3>{v.eye_position[0] - p[0], v.eye_position[1] - p[1], v.eye_position[2] - p[2]};
    };
    const std::array<double, 3> centre{0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
    for (std::size_t a = 0; a < 3; ++a) {
        // The hi face is turned toward the eye when the eye lies beyond its plane; away otherwise.
        std::array<double, 3> on_hi = centre;
        on_hi[a] = hi[a];
        const double d = toward_eye(on_hi)[a];
        out.away_hi[a] = d < -1e-9 * extent;  // edge-on: the lo side is the away one
    }
    const auto corner = [&](std::size_t a, std::size_t b, bool b_hi, std::size_t c, bool c_hi, bool at_hi) {
        std::array<double, 3> p{};
        p[a] = at_hi ? hi[a] : lo[a];
        p[b] = b_hi ? hi[b] : lo[b];
        p[c] = c_hi ? hi[c] : lo[c];
        return p;
    };
    const auto to_f3 = [](const std::array<double, 3>& p) {
        return Float3{static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])};
    };
    const D2D1_POINT_2F c_px = project(v, to_f3(centre), view);
    std::vector<std::array<float, 4>> all_kept;
    // The twelve edges: along axis a, at the sides (b_hi, c_hi) of the other two.
    for (std::size_t a = 0; a < 3; ++a) {
        const std::size_t b = (a + 1) % 3, c = (a + 2) % 3;
        struct Candidate {
            std::array<Float3, 2> edge;
            double nearness;
            float left;
        };
        std::vector<Candidate> outline;
        for (int k = 0; k < 4; ++k) {
            const bool b_hi = (k & 1) != 0, c_hi = (k & 2) != 0;
            const bool on_b = b_hi == out.away_hi[b], on_c = c_hi == out.away_hi[c];
            if (!on_b && !on_c) continue;  // on two faces toward the camera: not drawn
            const auto p0 = corner(a, b, b_hi, c, c_hi, false);
            const auto p1 = corner(a, b, b_hi, c, c_hi, true);
            const std::array<Float3, 2> e{to_f3(p0), to_f3(p1)};
            out.edges.push_back(e);
            if (on_b == on_c) continue;  // on two away faces: inside the frame
            const std::array<double, 3> mid{0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]), 0.5 * (p0[2] + p1[2])};
            const auto t = toward_eye(mid);
            const double nearness = v.perspective ? -std::sqrt(dotd(t, t)) : dotd(mid, t);
            outline.push_back({e, nearness, project(v, to_f3(mid), view).x});
        }
        if (outline.empty() || hi[a] - lo[a] <= 1e-9 * extent) continue;
        const Candidate& pick = *std::ranges::max_element(outline, [&](const Candidate& x, const Candidate& y) {
            return a == 1 ? x.left > y.left : x.nearness < y.nearness;  // y: leftmost; x, z: nearest
        });
        const D2D1_POINT_2F e0 = project(v, pick.edge[0], view), e1 = project(v, pick.edge[1], view);
        const float ex = e1.x - e0.x, ey = e1.y - e0.y;
        const float len = std::sqrt(ex * ex + ey * ey);
        if (len < 24.0f) continue;  // looking along the axis
        BoxAxis axis;
        axis.axis = static_cast<int>(a);
        axis.edge = pick.edge;
        // Outward: perpendicular to the edge on the screen, away from the box's centre.
        float ox = -ey / len, oy = ex / len;
        const float mx = 0.5f * (e0.x + e1.x) - c_px.x, my = 0.5f * (e0.y + e1.y) - c_px.y;
        if (ox * mx + oy * my < 0.0f) {
            ox = -ox;
            oy = -oy;
        }
        axis.outward = {ox, oy};
        std::vector<std::array<float, 4>> kept;  // this axis's kept labels: centres, half widths, half heights
        // Unit 26o (BSD3): a mapped axis ticks its own values, scene = origin + per_unit * value.
        const double origin = s.axes ? (*s.axes)[a].origin : 0.0;
        const double per_unit = s.axes ? (*s.axes)[a].per_unit : 1.0;
        const auto ticks = s.axes ? linear_ticks((lo[a] - origin) / per_unit, (hi[a] - origin) / per_unit, 5)
                                  : linear_ticks(lo[a], hi[a], 5);
        for (const Tick& t : ticks) {
            AxisLabel l;
            l.position = pick.edge[0];
            l.position[a] = static_cast<float>(s.axes ? origin + per_unit * t.value : t.value);
            const D2D1_POINT_2F p = project(v, l.position, view);
            l.x = p.x;
            l.y = p.y;
            l.value = t.value;
            l.text = t.label;
            const float half_w = 3.25f * static_cast<float>(t.label.size()), half_h = 7.0f;
            // Clear of the edge by the tick (6 px) and the label's own extent along the outward direction.
            const float clear = 10.0f + std::abs(ox) * half_w + std::abs(oy) * half_h;
            l.label_x = p.x + ox * clear;
            l.label_y = p.y + oy * clear;
            for (const auto& k : all_kept) {  // every axis's, so no two kept labels overlap
                const float dx = std::abs(l.label_x - k[0]), dy = std::abs(l.label_y - k[1]);
                if (dx < half_w + k[2] + 6.0f && dy < half_h + k[3] + 6.0f) l.kept = false;
            }
            if (l.kept) {
                kept.push_back({l.label_x, l.label_y, half_w, half_h});
                all_kept.push_back(kept.back());
            }
            axis.ticks.push_back(std::move(l));
        }
        static constexpr const char* titles[3] = {"x (\xC2\xB5m)", "y (\xC2\xB5m)", "z (\xC2\xB5m)"};
        axis.title = s.axes ? (*s.axes)[a].title : titles[a];
        float reach = 0.0f;
        for (const auto& k : kept) {
            reach = std::max(reach, (k[0] - 0.5f * (e0.x + e1.x)) * ox + (k[1] - 0.5f * (e0.y + e1.y)) * oy +
                                    std::abs(ox) * k[2] + std::abs(oy) * k[3]);
        }
        reach = std::max(reach, 16.0f);
        axis.title_x = 0.5f * (e0.x + e1.x) + ox * (reach + 14.0f);
        axis.title_y = 0.5f * (e0.y + e1.y) + oy * (reach + 14.0f);
        out.axes.push_back(std::move(axis));
    }
    return out;
}

// --- the view cube (Unit 26h, JD2, JD3) ------------------------------------------------------------

const std::vector<CameraPreset>& named_views() {
    static const std::vector<CameraPreset> views{
        {"Front", 0.0f, 0.0f},  {"Back", 180.0f, 0.0f}, {"Left", 270.0f, 0.0f},    {"Right", 90.0f, 0.0f},
        {"Top", 0.0f, 90.0f},   {"Bottom", 0.0f, -90.0f}, {"Isometric", 35.0f, 25.0f}};
    return views;
}

namespace {

const CameraPreset& named(std::string_view name) {
    for (const CameraPreset& p : named_views()) {
        if (p.name == name) return p;
    }
    return named_views().back();
}

// The isometric view from a cube corner: azimuth toward its x and z, elevation +-35.26 degrees.
CameraPreset corner_view(const std::array<int, 3>& c) {
    const double az = std::atan2(static_cast<double>(c[0]), static_cast<double>(c[2])) * 180.0 / std::numbers::pi;
    const double el = std::atan(1.0 / std::sqrt(2.0)) * 180.0 / std::numbers::pi;
    return {"Isometric", static_cast<float>(az), static_cast<float>(c[1] > 0 ? el : -el)};
}

bool inside_polygon(const std::array<std::array<float, 2>, 4>& p, float x, float y) {
    bool in = false;
    for (std::size_t i = 0, j = 3; i < 4; j = i++) {
        if ((p[i][1] > y) != (p[j][1] > y) &&
            x < (p[j][0] - p[i][0]) * (y - p[i][1]) / (p[j][1] - p[i][1]) + p[i][0]) {
            in = !in;
        }
    }
    return in;
}

}  // namespace

ViewCube view_cube(const FieldScene& s, int width, int height) {
    ViewCube out;
    if (s.dimension != 3 || s.look != Look::enhanced || !s.axis_triad) return out;
    const D2D1_RECT_F view = field_view_rect(s, width, height);
    const View v = make_view(s, (view.right - view.left) / (view.bottom - view.top));
    const std::array<double, 3> key = key_light(v);
    const std::array<double, 3> eye{-v.forward[0], -v.forward[1], -v.forward[2]};
    out.edge_px = 72.0f / std::sqrt(3.0f);
    out.centre = {view.left + 88.0f, view.bottom - 88.0f};  // room for an arrow and its label any way
    const float half = 0.5f * out.edge_px;
    const auto px = [&](const std::array<int, 3>& c) {
        const std::array<double, 3> p{double(c[0]), double(c[1]), double(c[2])};
        return std::array<float, 2>{out.centre[0] + half * static_cast<float>(dotd(p, v.right)),
                                    out.centre[1] - half * static_cast<float>(dotd(p, v.up))};
    };
    static constexpr const char* labels[3][2] = {{"Left", "Right"}, {"Bottom", "Top"}, {"Back", "Front"}};
    std::vector<std::array<int, 3>> shown_corners;
    for (int a = 0; a < 3; ++a) {
        for (const int sign : {1, -1}) {
            std::array<double, 3> n{};
            n[static_cast<std::size_t>(a)] = sign;
            if (dotd(n, eye) <= 1e-6) continue;  // turned away, or edge-on
            CubeFace f;
            f.label = labels[a][sign > 0 ? 1 : 0];
            f.axis = a;
            f.sign = sign;
            f.shade = 0.62f + 0.38f * static_cast<float>(std::max(0.0, dotd(n, key)));
            f.preset = named(f.label);
            // Around the face: the other two axes' signs in the order (-,-) (+,-) (+,+) (-,+).
            const int b = (a + 1) % 3, c = (a + 2) % 3;
            static constexpr int order[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
            for (int k = 0; k < 4; ++k) {
                std::array<int, 3> corner{};
                corner[static_cast<std::size_t>(a)] = sign;
                corner[static_cast<std::size_t>(b)] = order[k][0];
                corner[static_cast<std::size_t>(c)] = order[k][1];
                f.polygon[static_cast<std::size_t>(k)] = px(corner);
                if (std::ranges::find(shown_corners, corner) == shown_corners.end()) shown_corners.push_back(corner);
            }
            out.faces.push_back(std::move(f));
        }
    }
    for (const auto& c : shown_corners) out.corners.push_back({px(c), c, corner_view(c)});
    // The arrows from the shown corner nearest the screen's lower left, along the triad's directions.
    if (!out.corners.empty()) {
        const auto& origin = *std::ranges::min_element(out.corners, [](const CubeCorner& x, const CubeCorner& y) {
            return x.at[0] - x.at[1] < y.at[0] - y.at[1];
        });
        static constexpr const char* names[3] = {"x", "y", "z"};
        for (int a = 0; a < 3; ++a) {
            std::array<double, 3> e{};
            e[static_cast<std::size_t>(a)] = 1.0;
            const float dx = static_cast<float>(dotd(e, v.right)), dy = -static_cast<float>(dotd(e, v.up));
            const float len = out.edge_px + 10.0f;  // past the cube, even along its edges (head-on views)
            CubeArrow arrow{names[a], origin.at, origin.at};
            if (dx * dx + dy * dy > 1e-4f) arrow.to = {origin.at[0] + dx * len, origin.at[1] + dy * len};
            out.arrows.push_back(arrow);
        }
    }
    return out;
}

std::optional<CameraPreset> view_cube_hit(const FieldScene& s, int width, int height, float x, float y) {
    const ViewCube cube = view_cube(s, width, height);
    const CubeCorner* nearest = nullptr;
    float best = 7.0f * 7.0f;
    for (const CubeCorner& c : cube.corners) {
        const float dx = c.at[0] - x, dy = c.at[1] - y, d = dx * dx + dy * dy;
        if (d <= best) {
            best = d;
            nearest = &c;
        }
    }
    if (nearest) return nearest->preset;
    for (const CubeFace& f : cube.faces) {
        if (inside_polygon(f.polygon, x, y)) return f.preset;
    }
    return std::nullopt;
}

InsetFigure plane_inset_figure(const FieldScene& s, int width, int height) {
    InsetFigure out;
    if (s.dimension != 3 || !s.plane_inset) return out;
    const PlaneInset& p = *s.plane_inset;
    const D2D1_RECT_F view = field_view_rect(s, width, height);
    const View v = make_view(s, (view.right - view.left) / (view.bottom - view.top));
    constexpr float side = 140.0f, margin = 10.0f, pad = 12.0f;
    out.frame = {view.right - margin - side, view.top + margin, view.right - margin, view.top + margin + side};
    // The device's corners and the plane's, about the box's centre, on the camera's right and up.
    std::array<double, 3> c{};
    for (std::size_t k = 0; k < 3; ++k) c[k] = 0.5 * (p.lo_um[k] + p.hi_um[k]);
    const auto onto = [&](const std::array<double, 3>& q) {
        const std::array<double, 3> d{q[0] - c[0], q[1] - c[1], q[2] - c[2]};
        return std::array<double, 2>{dotd(d, v.right), dotd(d, v.up)};
    };
    std::array<std::array<double, 2>, 8> corner{};
    for (int i = 0; i < 8; ++i) {
        corner[static_cast<std::size_t>(i)] =
            onto({(i & 1) ? p.hi_um[0] : p.lo_um[0], (i & 2) ? p.hi_um[1] : p.lo_um[1], (i & 4) ? p.hi_um[2] : p.lo_um[2]});
    }
    const auto a = static_cast<std::size_t>(std::clamp(p.axis, 0, 2));
    const std::size_t b0 = (a + 1) % 3, b1 = (a + 2) % 3;
    std::array<std::array<double, 2>, 4> quad{};
    for (int i = 0; i < 4; ++i) {
        std::array<double, 3> q{};
        q[a] = std::clamp(p.at_um, p.lo_um[a], p.hi_um[a]);
        q[b0] = (i == 1 || i == 2) ? p.hi_um[b0] : p.lo_um[b0];
        q[b1] = (i >= 2) ? p.hi_um[b1] : p.lo_um[b1];
        quad[static_cast<std::size_t>(i)] = onto(q);
    }
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
    for (const auto& q : corner) {
        x0 = std::min(x0, q[0]);
        x1 = std::max(x1, q[0]);
        y0 = std::min(y0, q[1]);
        y1 = std::max(y1, q[1]);
    }
    const double span = std::max({x1 - x0, y1 - y0, 1e-12});
    const double k = (side - 2.0f * pad) / span;
    const float mx = 0.5f * (out.frame[0] + out.frame[2]), my = 0.5f * (out.frame[1] + out.frame[3]);
    const auto px = [&](const std::array<double, 2>& q) {
        return std::array<float, 2>{mx + static_cast<float>(k * (q[0] - 0.5 * (x0 + x1))),
                                    my - static_cast<float>(k * (q[1] - 0.5 * (y0 + y1)))};
    };
    for (int i = 0; i < 8; ++i) {
        for (int bit : {1, 2, 4}) {
            if (i & bit) continue;
            out.edges.push_back({px(corner[static_cast<std::size_t>(i)]), px(corner[static_cast<std::size_t>(i | bit)])});
        }
    }
    for (std::size_t i = 0; i < 4; ++i) out.plane[i] = px(quad[i]);
    return out;
}

}  // namespace NiTCAD::render
