#include "record.hpp"

#include <algorithm>
#include <bit>
#include <utility>

namespace NiTCAD::io::detail {

namespace {

void put(std::vector<std::byte>& out, std::uint64_t v, int bytes) {
    for (int b = 0; b < bytes; ++b) {
        out.push_back(static_cast<std::byte>((v >> (8 * b)) & 0xffU));
    }
}

std::uint64_t get(std::span<const std::byte> in, std::size_t at, int bytes) {
    std::uint64_t v = 0;
    for (int b = 0; b < bytes; ++b) {
        v |= static_cast<std::uint64_t>(std::to_integer<unsigned>(in[at + b])) << (8 * b);
    }
    return v;
}

}  // namespace

std::uint64_t fnv1a(std::span<const std::byte> bytes, std::uint64_t h) noexcept {
    for (const std::byte b : bytes) {
        h ^= std::to_integer<std::uint64_t>(b);
        h *= 0x100000001b3ULL;
    }
    return h;
}

base::Error corrupt(std::string what) {
    return {base::ErrorCode::invalid_input, "corrupt or incompatible file: " + std::move(what),
            std::nullopt};
}

void Writer::field(std::uint16_t key, std::span<const std::byte> payload) {
    put(bytes_, key, 2);
    put(bytes_, payload.size(), 8);
    bytes_.insert(bytes_.end(), payload.begin(), payload.end());
}

void Writer::u64(std::uint16_t key, std::uint64_t v) {
    std::vector<std::byte> b;
    put(b, v, 8);
    field(key, b);
}

void Writer::i64(std::uint16_t key, std::int64_t v) { u64(key, static_cast<std::uint64_t>(v)); }

void Writer::f64(std::uint16_t key, double v) { u64(key, std::bit_cast<std::uint64_t>(v)); }

void Writer::text(std::uint16_t key, std::string_view v) {
    field(key, std::as_bytes(std::span(v.data(), v.size())));
}

void Writer::f64s(std::uint16_t key, std::span<const double> v) {
    std::vector<std::byte> b;
    b.reserve(8 * v.size());
    for (const double d : v) put(b, std::bit_cast<std::uint64_t>(d), 8);
    field(key, b);
}

void Writer::i64s(std::uint16_t key, std::span<const std::int64_t> v) {
    std::vector<std::byte> b;
    b.reserve(8 * v.size());
    for (const std::int64_t d : v) put(b, static_cast<std::uint64_t>(d), 8);
    field(key, b);
}

void Writer::record(std::uint16_t key, const Writer& nested) { field(key, nested.bytes_); }

std::expected<Reader, base::Error> Reader::parse(std::span<const std::byte> bytes) {
    Reader r;
    std::size_t at = 0;
    while (at < bytes.size()) {
        if (bytes.size() - at < 10) return std::unexpected(corrupt("truncated field header"));
        const auto key = static_cast<std::uint16_t>(get(bytes, at, 2));
        const std::uint64_t length = get(bytes, at + 2, 8);
        at += 10;
        if (length > bytes.size() - at) return std::unexpected(corrupt("field past the end"));
        r.fields_.push_back({key, bytes.subspan(at, static_cast<std::size_t>(length))});
        at += static_cast<std::size_t>(length);
    }
    return r;
}

std::optional<std::span<const std::byte>> Reader::find(std::uint16_t key) const noexcept {
    for (const Field& f : fields_) {
        if (f.key == key) return f.bytes;
    }
    return std::nullopt;
}

bool Reader::has(std::uint16_t key) const noexcept { return find(key).has_value(); }

std::expected<std::span<const std::byte>, base::Error> Reader::need(std::uint16_t key) const {
    if (auto b = find(key)) return *b;
    return std::unexpected(corrupt("missing field " + std::to_string(key)));
}

std::expected<std::uint64_t, base::Error> Reader::u64(std::uint16_t key) const {
    auto b = need(key);
    if (!b) return std::unexpected(std::move(b.error()));
    if (b->size() != 8) return std::unexpected(corrupt("field " + std::to_string(key) + " size"));
    return get(*b, 0, 8);
}

std::expected<std::int64_t, base::Error> Reader::i64(std::uint16_t key) const {
    auto v = u64(key);
    if (!v) return std::unexpected(std::move(v.error()));
    return static_cast<std::int64_t>(*v);
}

std::expected<double, base::Error> Reader::f64(std::uint16_t key) const {
    auto v = u64(key);
    if (!v) return std::unexpected(std::move(v.error()));
    return std::bit_cast<double>(*v);
}

std::expected<std::string, base::Error> Reader::text(std::uint16_t key) const {
    auto b = need(key);
    if (!b) return std::unexpected(std::move(b.error()));
    std::string s(b->size(), '\0');
    std::transform(b->begin(), b->end(), s.begin(),
                   [](std::byte c) { return static_cast<char>(std::to_integer<unsigned>(c)); });
    return s;
}

std::expected<std::vector<double>, base::Error> Reader::f64s(std::uint16_t key) const {
    auto b = need(key);
    if (!b) return std::unexpected(std::move(b.error()));
    if (b->size() % 8 != 0) {
        return std::unexpected(corrupt("field " + std::to_string(key) + " size"));
    }
    std::vector<double> v(b->size() / 8);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = std::bit_cast<double>(get(*b, 8 * i, 8));
    return v;
}

std::expected<std::vector<std::int64_t>, base::Error> Reader::i64s(std::uint16_t key) const {
    auto b = need(key);
    if (!b) return std::unexpected(std::move(b.error()));
    if (b->size() % 8 != 0) {
        return std::unexpected(corrupt("field " + std::to_string(key) + " size"));
    }
    std::vector<std::int64_t> v(b->size() / 8);
    for (std::size_t i = 0; i < v.size(); ++i) {
        v[i] = static_cast<std::int64_t>(get(*b, 8 * i, 8));
    }
    return v;
}

std::expected<Reader, base::Error> Reader::record(std::uint16_t key) const {
    auto b = need(key);
    if (!b) return std::unexpected(std::move(b.error()));
    return parse(*b);
}

std::expected<std::vector<Reader>, base::Error> Reader::records(std::uint16_t key) const {
    std::vector<Reader> out;
    for (const Field& f : fields_) {
        if (f.key != key) continue;
        auto r = parse(f.bytes);
        if (!r) return std::unexpected(std::move(r.error()));
        out.push_back(std::move(*r));
    }
    return out;
}

std::expected<std::uint64_t, base::Error> Reader::u64_or(std::uint16_t key,
                                                        std::uint64_t fallback) const {
    return has(key) ? u64(key) : fallback;
}

std::expected<std::int64_t, base::Error> Reader::i64_or(std::uint16_t key,
                                                       std::int64_t fallback) const {
    return has(key) ? i64(key) : fallback;
}

std::expected<std::vector<double>, base::Error> Reader::f64s_or_empty(std::uint16_t key) const {
    return has(key) ? f64s(key) : std::vector<double>{};
}

}  // namespace NiTCAD::io::detail
