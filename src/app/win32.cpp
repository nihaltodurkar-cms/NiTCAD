// Win32 helpers (win32.hpp).
#include "win32.hpp"

#include <commdlg.h>
#include <dxgi1_3.h>
#include <dxgidebug.h>
#include <shlobj.h>
#include <wrl/client.h>

#include <vector>

namespace NiTCAD::app::win32 {

std::wstring widen(std::string_view s) {
    if (s.empty()) return {};
    const int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
    std::wstring w(static_cast<std::size_t>(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), w.data(), n);
    return w;
}

std::string narrow(std::wstring_view w) {
    if (w.empty()) return {};
    const int n = WideCharToMultiByte(CP_UTF8, 0, w.data(), static_cast<int>(w.size()), nullptr, 0, nullptr, nullptr);
    std::string s(static_cast<std::size_t>(n), '\0');
    WideCharToMultiByte(CP_UTF8, 0, w.data(), static_cast<int>(w.size()), s.data(), n, nullptr, nullptr);
    return s;
}

std::string text_of(HWND window) {
    const int n = GetWindowTextLengthW(window);
    std::wstring w(static_cast<std::size_t>(n) + 1, L'\0');
    const int got = GetWindowTextW(window, w.data(), n + 1);
    w.resize(static_cast<std::size_t>(std::max(got, 0)));
    return narrow(w);
}

void set_text(HWND window, std::string_view utf8) { SetWindowTextW(window, widen(utf8).c_str()); }

HFONT ui_font(int dpi, bool bold) {
    return CreateFontW(-MulDiv(9, dpi, 72), 0, 0, 0, bold ? FW_SEMIBOLD : FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
                       OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
}

HFONT code_font(int dpi) {
    return CreateFontW(-MulDiv(10, dpi, 72), 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                       CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, FIXED_PITCH | FF_MODERN, L"Consolas");
}

int dpi_of(HWND window) {
    const UINT d = window ? GetDpiForWindow(window) : 0;
    return d ? static_cast<int>(d) : 96;
}

bool register_class(const wchar_t* name, WNDPROC proc, HBRUSH background, HCURSOR cursor, UINT style) {
    WNDCLASSEXW existing{sizeof(existing)};
    if (GetClassInfoExW(GetModuleHandleW(nullptr), name, &existing)) return true;
    WNDCLASSEXW c{sizeof(c)};
    c.style = style;
    c.lpfnWndProc = proc;
    c.hInstance = GetModuleHandleW(nullptr);
    c.hCursor = cursor ? cursor : LoadCursorW(nullptr, IDC_ARROW);
    c.hbrBackground = background;
    c.lpszClassName = name;
    c.hIcon = LoadIconW(nullptr, IDI_APPLICATION);
    return RegisterClassExW(&c) != 0;
}

std::filesystem::path state_file() {
    PWSTR folder = nullptr;
    std::filesystem::path out;
    if (SUCCEEDED(SHGetKnownFolderPath(FOLDERID_RoamingAppData, 0, nullptr, &folder))) {
        out = std::filesystem::path(folder) / L"NiTCAD" / L"app.ini";
    }
    CoTaskMemFree(folder);
    return out;
}

std::filesystem::path documents_dir() {
    PWSTR folder = nullptr;
    std::filesystem::path out;
    if (SUCCEEDED(SHGetKnownFolderPath(FOLDERID_Documents, 0, nullptr, &folder))) out = folder;
    CoTaskMemFree(folder);
    return out;
}

int live_graphics_objects() {
    Microsoft::WRL::ComPtr<IDXGIDebug1> debug;
    Microsoft::WRL::ComPtr<IDXGIInfoQueue> queue;
    if (FAILED(DXGIGetDebugInterface1(0, IID_PPV_ARGS(&debug))) || FAILED(DXGIGetDebugInterface1(0, IID_PPV_ARGS(&queue)))) {
        return -1;
    }
    queue->ClearStoredMessages(DXGI_DEBUG_ALL);
    debug->ReportLiveObjects(DXGI_DEBUG_ALL,
                             static_cast<DXGI_DEBUG_RLO_FLAGS>(DXGI_DEBUG_RLO_DETAIL | DXGI_DEBUG_RLO_IGNORE_INTERNAL));
    int live = 0;
    const UINT64 n = queue->GetNumStoredMessages(DXGI_DEBUG_ALL);
    for (UINT64 i = 0; i < n; ++i) {
        SIZE_T size = 0;
        queue->GetMessage(DXGI_DEBUG_ALL, i, nullptr, &size);
        std::vector<char> buffer(size);
        auto* m = reinterpret_cast<DXGI_INFO_QUEUE_MESSAGE*>(buffer.data());
        if (SUCCEEDED(queue->GetMessage(DXGI_DEBUG_ALL, i, m, &size)) && std::string_view(m->pDescription).starts_with("Live ")) {
            ++live;
        }
    }
    return live;
}

std::filesystem::path open_dialog(HWND owner, const wchar_t* filter, const wchar_t* title) {
    std::vector<wchar_t> buffer(32768, L'\0');
    OPENFILENAMEW o{sizeof(o)};
    o.hwndOwner = owner;
    o.lpstrFilter = filter;
    o.lpstrFile = buffer.data();
    o.nMaxFile = static_cast<DWORD>(buffer.size());
    o.lpstrTitle = title;
    o.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_EXPLORER;
    if (!GetOpenFileNameW(&o)) return {};
    return std::filesystem::path(buffer.data());
}

std::filesystem::path save_dialog(HWND owner, const wchar_t* filter, const wchar_t* title, const wchar_t* extension) {
    std::vector<wchar_t> buffer(32768, L'\0');
    OPENFILENAMEW o{sizeof(o)};
    o.hwndOwner = owner;
    o.lpstrFilter = filter;
    o.lpstrFile = buffer.data();
    o.nMaxFile = static_cast<DWORD>(buffer.size());
    o.lpstrTitle = title;
    o.lpstrDefExt = extension;
    o.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST | OFN_EXPLORER;
    if (!GetSaveFileNameW(&o)) return {};
    return std::filesystem::path(buffer.data());
}

}  // namespace NiTCAD::app::win32
