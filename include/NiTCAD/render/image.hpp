// Drawing scenes into images (Unit 25a; DECISIONS.md VR5, VR8, VQ1, VQ2, VQ4, VG4, VG7), and into a
// window's swap chain (Unit 26, VA4).
//
// A Renderer owns a Direct3D 12 device (a hardware adapter, or WARP, the software rasterizer that
// ships with Windows) and draws a FieldScene or a PlotScene off screen:
// - triangles and lines on Direct3D 12, at twice the image's resolution and averaged down (2 x 2
//   supersampling); field triangles look the colormap up per pixel from the interpolated position;
//   lit triangles take one directional light and ambient, two-sided; translucent triangles go through
//   weighted blended order-independent transparency (McGuire and Bavoil 2013), so their order does
//   not matter; lines are screen-space quads of their width;
// - text, the colorbar, axes, the axis triad and plots on Direct2D and DirectWrite, drawn into the
//   same Direct3D 12 target through D3D11On12 (VQ1).
// Shaders are HLSL compiled at build time with the SDK's dxc and embedded (VQ2).
//
// The camera of a FieldScene sets the direction and the projection; drawing fits the scene's bounds
// to the image's aspect. A Renderer is not safe to use from two threads at once.
//
// Errors: resource_exhausted when the graphics device cannot be created or fails (the message
// names the call and its HRESULT; the context value is the HRESULT), or memory runs out;
// invalid_input for an image size outside 16 to 16384 pixels.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/render/color.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::render {

// 8-bit sRGB with alpha, rows top to bottom, 4 bytes per pixel (r, g, b, a).
struct Image {
    int width = 0;
    int height = 0;
    std::vector<std::uint8_t> rgba;
};

enum class Adapter : std::uint8_t {
    automatic,  // the first hardware adapter that supports Direct3D 12, else WARP
    hardware,
    warp,
};

struct RendererOptions {
    Adapter adapter = Adapter::automatic;
    // The Direct3D 12 debug layer (needs the Graphics Tools feature); its error messages are
    // counted by debug_errors().
    bool debug_layer = false;
};

// A window's swap chain on a renderer's queue (Unit 26, DECISIONS.md VA4): Renderer::present draws a
// scene and shows it in the window, the very image draw() would return at the window's size (the
// image is copied into the swap chain's back buffer, so the window and a PNG of the view agree).
// Made by Renderer::window_target for a window handle (an HWND); it must not outlive its renderer.
class WindowTarget {
public:
    WindowTarget(WindowTarget&&) noexcept;
    WindowTarget& operator=(WindowTarget&&) noexcept;
    WindowTarget(const WindowTarget&) = delete;
    WindowTarget& operator=(const WindowTarget&) = delete;
    ~WindowTarget();

    [[nodiscard]] int width() const noexcept;
    [[nodiscard]] int height() const noexcept;

    struct Impl;  // the swap chain (renderer.cpp)

private:
    friend class Renderer;
    explicit WindowTarget(std::unique_ptr<Impl> impl) noexcept;
    std::unique_ptr<Impl> impl_;
};

// Where a world point (um) of a FieldScene lands in an image of the scene drawn at this size: the
// pixel (x right, y down) and the depth in [0, 1] (nearer smaller), as draw() lays the view out.
struct ScreenPoint {
    float x = 0.0f;
    float y = 0.0f;
    float depth = 0.0f;
};
[[nodiscard]] ScreenPoint screen_position(const FieldScene& scene, int width, int height, const Float3& point);

class Renderer {
public:
    [[nodiscard]] static std::expected<Renderer, base::Error> create(const RendererOptions& options = {});

    // A swap chain for a window (an HWND) of the given client size (16 to 16384 pixels). Errors:
    // resource_exhausted when DXGI cannot make it (no desktop, a window that is not a window).
    [[nodiscard]] std::expected<WindowTarget, base::Error> window_target(void* window, int width, int height);
    // Draws the scene at the size (resizing the swap chain when it changed) and presents it.
    [[nodiscard]] std::expected<void, base::Error> present(const FieldScene& scene, WindowTarget& target, int width,
                                                           int height);
    [[nodiscard]] std::expected<void, base::Error> present(const PlotScene& scene, WindowTarget& target, int width,
                                                           int height);

    Renderer(Renderer&&) noexcept;
    Renderer& operator=(Renderer&&) noexcept;
    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;
    ~Renderer();

    [[nodiscard]] std::expected<Image, base::Error> draw(const FieldScene& scene, int width, int height);
    [[nodiscard]] std::expected<Image, base::Error> draw(const PlotScene& scene, int width, int height);

    // The time per frame of the scene's Direct3D passes (opaque, lines, transparency, the 2 x 2
    // resolve), repeated `frames` times on geometry already on the GPU, in milliseconds: what a
    // window redraws cost (Unit 26), without the per-image upload, overlays and readback of draw.
    [[nodiscard]] std::expected<double, base::Error> frame_time(const FieldScene& scene, int width,
                                                               int height, int frames);

    // The adapter's description, and whether it is WARP.
    [[nodiscard]] const std::string& adapter_name() const noexcept;
    [[nodiscard]] bool is_warp() const noexcept;
    // Error and corruption messages the debug layer has stored so far (0 without the layer).
    [[nodiscard]] std::size_t debug_errors() const;
    // Their descriptions.
    [[nodiscard]] std::vector<std::string> debug_messages() const;

private:
    struct Impl;
    explicit Renderer(std::unique_ptr<Impl> impl) noexcept;
    std::unique_ptr<Impl> impl_;
};

// PNG files through WIC. Errors: resource_exhausted when the file cannot be written or read (the
// message names the path), invalid_input for an image whose size and data disagree.
[[nodiscard]] std::expected<void, base::Error> save_png(const Image& image,
                                                       const std::filesystem::path& path);
[[nodiscard]] std::expected<Image, base::Error> load_png(const std::filesystem::path& path);

// How two images of the same size differ: the largest difference of any channel (0 to 255) and the
// fraction of pixels with a channel differing by more than `tolerance`.
struct ImageDifference {
    int max_channel = 0;
    double fraction_over = 0.0;
};
[[nodiscard]] std::expected<ImageDifference, base::Error> compare(const Image& a, const Image& b,
                                                                 int tolerance);

// Images of one size tiled in rows of `columns`, left to right and top to bottom, on white (Unit 25c,
// DC4: a filmstrip). Errors (invalid_input): no image, images of different sizes, columns below 1.
[[nodiscard]] std::expected<Image, base::Error> tile_images(std::span<const Image> images, int columns);

}  // namespace NiTCAD::render
