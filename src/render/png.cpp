// PNG files through WIC, and image comparison (image.hpp).
#include "NiTCAD/render/image.hpp"

#include <windows.h>
#include <wincodec.h>
#include <wrl/client.h>

#include <algorithm>
#include <cstdlib>
#include <format>
#include <string>

namespace NiTCAD::render {

namespace {

using Microsoft::WRL::ComPtr;

base::Error file_error(std::string_view what, const std::filesystem::path& path, HRESULT hr) {
    return {base::ErrorCode::resource_exhausted,
            std::format("{} '{}' failed (HRESULT 0x{:08X})", what, path.string(), static_cast<unsigned>(hr)),
            base::ErrorContext{.index = std::nullopt, .value = static_cast<double>(static_cast<unsigned>(hr))}};
}

// COM for this thread, initialized once (apartment-threaded, as WIC allows).
bool com_ready() {
    thread_local const HRESULT hr = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
    return SUCCEEDED(hr) || hr == RPC_E_CHANGED_MODE;
}

std::expected<ComPtr<IWICImagingFactory>, HRESULT> wic() {
    if (!com_ready()) return std::unexpected(E_FAIL);
    ComPtr<IWICImagingFactory> f;
    const HRESULT hr = CoCreateInstance(CLSID_WICImagingFactory, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&f));
    if (FAILED(hr)) return std::unexpected(hr);
    return f;
}

}  // namespace

std::expected<void, base::Error> save_png(const Image& image, const std::filesystem::path& path) {
    if (image.width <= 0 || image.height <= 0 ||
        image.rgba.size() != static_cast<std::size_t>(image.width) * image.height * 4) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           "the image's size and data disagree", std::nullopt});
    }
    auto f = wic();
    if (!f) return std::unexpected(file_error("creating WIC for", path, f.error()));
    ComPtr<IWICStream> stream;
    HRESULT hr = (*f)->CreateStream(&stream);
    if (SUCCEEDED(hr)) hr = stream->InitializeFromFilename(path.c_str(), GENERIC_WRITE);
    if (FAILED(hr)) return std::unexpected(file_error("opening", path, hr));
    ComPtr<IWICBitmapEncoder> encoder;
    hr = (*f)->CreateEncoder(GUID_ContainerFormatPng, nullptr, &encoder);
    if (SUCCEEDED(hr)) hr = encoder->Initialize(stream.Get(), WICBitmapEncoderNoCache);
    ComPtr<IWICBitmapFrameEncode> frame;
    if (SUCCEEDED(hr)) hr = encoder->CreateNewFrame(&frame, nullptr);
    if (SUCCEEDED(hr)) hr = frame->Initialize(nullptr);
    if (SUCCEEDED(hr)) hr = frame->SetSize(static_cast<UINT>(image.width), static_cast<UINT>(image.height));
    // The PNG encoder takes BGRA, not RGBA.
    WICPixelFormatGUID format = GUID_WICPixelFormat32bppBGRA;
    if (SUCCEEDED(hr)) hr = frame->SetPixelFormat(&format);
    if (SUCCEEDED(hr) && format != GUID_WICPixelFormat32bppBGRA) hr = WINCODEC_ERR_UNSUPPORTEDPIXELFORMAT;
    if (SUCCEEDED(hr)) {
        std::vector<std::uint8_t> bgra(image.rgba);
        for (std::size_t p = 0; p + 3 < bgra.size(); p += 4) std::swap(bgra[p], bgra[p + 2]);
        hr = frame->WritePixels(static_cast<UINT>(image.height), static_cast<UINT>(image.width) * 4,
                                static_cast<UINT>(bgra.size()), bgra.data());
    }
    if (SUCCEEDED(hr)) hr = frame->Commit();
    if (SUCCEEDED(hr)) hr = encoder->Commit();
    if (FAILED(hr)) return std::unexpected(file_error("writing PNG", path, hr));
    return {};
}

std::expected<Image, base::Error> load_png(const std::filesystem::path& path) {
    auto f = wic();
    if (!f) return std::unexpected(file_error("creating WIC for", path, f.error()));
    ComPtr<IWICBitmapDecoder> decoder;
    HRESULT hr = (*f)->CreateDecoderFromFilename(path.c_str(), nullptr, GENERIC_READ,
                                                 WICDecodeMetadataCacheOnDemand, &decoder);
    ComPtr<IWICBitmapFrameDecode> frame;
    if (SUCCEEDED(hr)) hr = decoder->GetFrame(0, &frame);
    ComPtr<IWICFormatConverter> converter;
    if (SUCCEEDED(hr)) hr = (*f)->CreateFormatConverter(&converter);
    if (SUCCEEDED(hr)) {
        hr = converter->Initialize(frame.Get(), GUID_WICPixelFormat32bppRGBA, WICBitmapDitherTypeNone, nullptr, 0.0,
                                   WICBitmapPaletteTypeCustom);
    }
    UINT w = 0, h = 0;
    if (SUCCEEDED(hr)) hr = converter->GetSize(&w, &h);
    Image img{static_cast<int>(w), static_cast<int>(h), std::vector<std::uint8_t>(static_cast<std::size_t>(w) * h * 4)};
    if (SUCCEEDED(hr)) hr = converter->CopyPixels(nullptr, w * 4, static_cast<UINT>(img.rgba.size()), img.rgba.data());
    if (FAILED(hr)) return std::unexpected(file_error("reading PNG", path, hr));
    return img;
}

std::expected<ImageDifference, base::Error> compare(const Image& a, const Image& b, int tolerance) {
    if (a.width != b.width || a.height != b.height || a.rgba.size() != b.rgba.size()) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, "the images differ in size", std::nullopt});
    }
    ImageDifference d;
    std::size_t over = 0;
    for (std::size_t p = 0; p + 3 < a.rgba.size(); p += 4) {
        int worst = 0;
        for (std::size_t c = 0; c < 4; ++c) worst = std::max(worst, std::abs(int{a.rgba[p + c]} - int{b.rgba[p + c]}));
        d.max_channel = std::max(d.max_channel, worst);
        if (worst > tolerance) ++over;
    }
    const std::size_t pixels = a.rgba.size() / 4;
    d.fraction_over = pixels ? static_cast<double>(over) / static_cast<double>(pixels) : 0.0;
    return d;
}

std::expected<Image, base::Error> tile_images(std::span<const Image> images, int columns) {
    const auto invalid = [](const char* m) { return base::Error{base::ErrorCode::invalid_input, m, std::nullopt}; };
    if (images.empty()) return std::unexpected(invalid("no image to tile"));
    if (columns < 1) return std::unexpected(invalid("tiles need at least one column"));
    const int w = images[0].width, h = images[0].height;
    for (const Image& i : images) {
        if (i.width != w || i.height != h || i.rgba.size() != static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 4) {
            return std::unexpected(invalid("the images to tile differ in size"));
        }
    }
    const int cols = std::min(columns, static_cast<int>(images.size()));
    const int rows = (static_cast<int>(images.size()) + cols - 1) / cols;
    Image out;
    out.width = w * cols;
    out.height = h * rows;
    out.rgba.assign(static_cast<std::size_t>(out.width) * static_cast<std::size_t>(out.height) * 4, 255);
    for (std::size_t k = 0; k < images.size(); ++k) {
        const int x0 = static_cast<int>(k) % cols * w, y0 = static_cast<int>(k) / cols * h;
        for (int y = 0; y < h; ++y) {
            const auto src = images[k].rgba.begin() + static_cast<std::ptrdiff_t>(y) * w * 4;
            std::copy(src, src + static_cast<std::ptrdiff_t>(w) * 4,
                      out.rgba.begin() + (static_cast<std::ptrdiff_t>(y0 + y) * out.width + x0) * 4);
        }
    }
    return out;
}

}  // namespace NiTCAD::render
