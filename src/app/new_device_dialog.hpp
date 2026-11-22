// File > New device... (Unit 26n slice 1, DT4, DECISIONS.md): an owned popup over the main window, which is
// disabled while it is open (no nested message loop). The library's families in a list; the chosen deck's
// thumbnail and its main dimensions (the first [parameters] group, a DimensionsPanel over the deck's text in
// memory); the folder and the name; Create writes the deck there (never over a file) and hands its path to the
// host, which opens it on the Geometry tab. The model is app/new_device.hpp.
#pragma once

#include "NiTCAD/app/new_device.hpp"
#include "dimensions_panel.hpp"

#include <windows.h>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace NiTCAD::app {

class NewDeviceDialog {
public:
    struct Host {
        virtual void new_device_created(const std::filesystem::path& deck) = 0;  // the host closes the dialog later (posted)
        virtual void new_device_closed() = 0;                                    // Cancel or the close box

    protected:
        ~Host() = default;
    };

    NewDeviceDialog(HWND owner, Host& host, std::vector<devices::Family> families, std::filesystem::path folder, int dpi,
                    HFONT font, HFONT bold);
    ~NewDeviceDialog();
    NewDeviceDialog(const NewDeviceDialog&) = delete;
    NewDeviceDialog& operator=(const NewDeviceDialog&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return hwnd_; }
    [[nodiscard]] const std::vector<devices::Family>& families() const noexcept { return families_; }
    void select(std::size_t family);  // its deck's text, its name, its main dimensions
    [[nodiscard]] std::size_t selected() const noexcept { return selected_; }
    [[nodiscard]] DimensionsPanel& panel() noexcept { return *panel_; }
    [[nodiscard]] const std::string& text() const noexcept { return text_; }
    void set_folder(const std::filesystem::path& folder);
    [[nodiscard]] std::filesystem::path folder() const;
    void set_name(const std::string& name);
    [[nodiscard]] std::string name() const;
    [[nodiscard]] bool can_create() const noexcept { return thumb_.error.empty() && !text_.empty(); }
    [[nodiscard]] const std::string& message() const noexcept { return message_; }
    bool create();  // as the Create button
    void cancel();  // as Cancel

private:
    struct PanelHost final : DimensionsPanel::Host {
        NewDeviceDialog* d = nullptr;
        std::string dimensions_deck() override { return d->text_; }
        void dimensions_edited(const std::string& text) override;
        void dimensions_selected(const std::string&) override {}
        std::string dimensions_refuse(const std::string& text) override;
    };
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK thumb_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    void layout();
    void show_message();
    void paint_thumbnail(HDC dc, const RECT& r) const;

    HWND owner_ = nullptr;
    Host& host_;
    HWND hwnd_ = nullptr, list_ = nullptr, thumb_wnd_ = nullptr, message_label_ = nullptr;
    HWND folder_label_ = nullptr, folder_ = nullptr, browse_ = nullptr, name_label_ = nullptr, name_ = nullptr;
    HWND create_ = nullptr, cancel_ = nullptr;
    PanelHost panel_host_;
    std::unique_ptr<DimensionsPanel> panel_;
    std::vector<devices::Family> families_;
    std::size_t selected_ = 0;
    std::string text_, message_;
    devices::Thumbnail thumb_;
    int dpi_ = 96;
    HFONT font_ = nullptr, bold_ = nullptr;
};

}  // namespace NiTCAD::app
