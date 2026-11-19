// The device and curve view of the main window (Unit 26, DECISIONS.md VA4, VA8): a child window that
// presents a render scene in a DXGI swap chain on the renderer's queue, and turns the mouse and keys into
// camera moves, probes and fits for its host.
//
// Mouse: left drag orbits (3D) or pans (2D), right or middle drag pans, the wheel zooms, a double click
// fits, a click (no drag) probes. Keys while focused: arrows orbit (2D: pan), Shift+arrows pan, + and -
// zoom, Home fits.
#pragma once

#include <windows.h>

#include <cstddef>
#include <optional>
#include <string>

#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/render/image.hpp"

namespace NiTCAD::app {

class Viewport {
public:
    struct Host {
        virtual void viewport_orbit(float azimuth_deg, float elevation_deg) = 0;
        virtual void viewport_pan(float dx, float dy) = 0;  // fractions of the half view
        virtual void viewport_zoom(float factor) = 0;
        virtual void viewport_fit() = 0;
        virtual void viewport_click(float x, float y) = 0;  // client pixels
        virtual void viewport_hover(float x, float y) = 0;
        // The mouse left the view, or a drag began (Unit 26e, HD3).
        virtual void viewport_leave() = 0;
        virtual bool viewport_is_3d() const = 0;
        // Whether a click there is the view cube's (Unit 26h, JD3): the cursor becomes a hand.
        virtual bool viewport_over_cube(float x, float y) const = 0;
        // A numpad key's named view (Unit 26h, JD4): an index into render::named_views().
        virtual void viewport_named_view(int index) = 0;

    protected:
        ~Host() = default;
    };

    Viewport(HWND parent, Host& host, int id);
    ~Viewport();
    Viewport(const Viewport&) = delete;
    Viewport& operator=(const Viewport&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return hwnd_; }
    // What to show: a scene drawn with the renderer (both outlive the next call), or a message.
    void show(render::Renderer* renderer, const Scene* scene);
    void show_message(std::string message);
    // Makes the swap chain now (at the view's size), so a first present later does not pay for it: the
    // run monitor prepares its plot before a run starts.
    void prepare(render::Renderer* renderer);
    // Draws now (also on WM_PAINT and WM_SIZE). Returns the error of the last present, if any.
    void present();
    // Unit 26i (ND2): draws again after a change of the scene's camera only, reusing the swap chain's
    // buffers (Renderer::present_camera). A plot, or a first frame, is drawn as present() draws it.
    void present_camera();
    // Draws a pending camera frame now (tests and the benchmark; the window draws it on WM_PAINT).
    void flush();
    // Frames presented so far (NG3).
    [[nodiscard]] std::size_t presents() const noexcept { return presents_; }
    [[nodiscard]] const std::optional<std::string>& error() const noexcept { return error_; }
    [[nodiscard]] double last_present_ms() const noexcept { return last_ms_; }
    [[nodiscard]] int width() const noexcept;
    [[nodiscard]] int height() const noexcept;
    // The swap chain is released (before its renderer is).
    void release();

private:
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    LRESULT handle(UINT message, WPARAM w, LPARAM l);
    void paint_message(HDC dc);
    void draw(bool camera_only);

    HWND hwnd_ = nullptr;
    Host& host_;
    render::Renderer* renderer_ = nullptr;
    const Scene* scene_ = nullptr;
    std::optional<render::WindowTarget> target_;
    std::string message_;
    std::optional<std::string> error_;
    double last_ms_ = 0.0;
    // Dragging.
    bool dragging_ = false;
    bool moved_ = false;
    bool tracking_ = false;  // TrackMouseEvent's leave is pending
    bool camera_pending_ = false;  // a camera frame waits for WM_PAINT (Unit 26i, ND4)
    std::size_t presents_ = 0;
    int button_ = 0;  // 1 left, 2 right, 3 middle
    POINT last_{};
    POINT down_{};
};

}  // namespace NiTCAD::app
