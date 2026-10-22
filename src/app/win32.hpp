// Small Win32 helpers of the application's window code (Unit 26): text between UTF-8 and UTF-16, fonts
// at a monitor's DPI, and window classes.
#pragma once

#include <windows.h>

#include <filesystem>
#include <string>
#include <string_view>

namespace NiTCAD::app::win32 {

[[nodiscard]] std::wstring widen(std::string_view utf8);
[[nodiscard]] std::string narrow(std::wstring_view utf16);

// A window's text as UTF-8, and setting it.
[[nodiscard]] std::string text_of(HWND window);
void set_text(HWND window, std::string_view utf8);

// The UI font and the deck's fixed-pitch font at a DPI (points: 9 and 10). The caller deletes them.
[[nodiscard]] HFONT ui_font(int dpi, bool bold = false);
[[nodiscard]] HFONT code_font(int dpi);

// A window's DPI (96 when unknown).
[[nodiscard]] int dpi_of(HWND window);

// Registers a window class once per process (by name); returns false when that fails.
bool register_class(const wchar_t* name, WNDPROC proc, HBRUSH background, HCURSOR cursor = nullptr,
                    UINT style = CS_HREDRAW | CS_VREDRAW);

// The user's application-data file for the state (%APPDATA%\NiTCAD\app.ini), empty when unknown.
[[nodiscard]] std::filesystem::path state_file();

// Open and save dialogs: the chosen path, or empty when cancelled. `filter` is a double-NUL list.
[[nodiscard]] std::filesystem::path open_dialog(HWND owner, const wchar_t* filter, const wchar_t* title);
[[nodiscard]] std::filesystem::path save_dialog(HWND owner, const wchar_t* filter, const wchar_t* title,
                                                const wchar_t* default_extension);

}  // namespace NiTCAD::app::win32
