// The record encoding of the io layer (private; ARCHITECTURE.md 6.5, 6.6 and 11.3, Unit 29).
//
// A record is a sequence of fields, each a 16-bit key, a 64-bit byte length and the bytes; a key
// may repeat (a list is its key repeated, in order). Every integer is little-endian whatever the
// host, and a double is its IEEE-754 bit pattern, so a value read back is the value written, bit
// for bit. A reader looks fields up by key: it skips keys it does not know (a newer writer's) and
// takes a default for keys that are absent (an older writer's), so the format can grow without
// breaking either side. Nested records are fields whose bytes are a record.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::io::detail {

// 64-bit FNV-1a over bytes, continuing from `h`.
std::uint64_t fnv1a(std::span<const std::byte> bytes,
                    std::uint64_t h = 0xcbf29ce484222325ULL) noexcept;

class Writer {
public:
    void u64(std::uint16_t key, std::uint64_t v);
    void i64(std::uint16_t key, std::int64_t v);
    void f64(std::uint16_t key, double v);
    void text(std::uint16_t key, std::string_view v);
    void f64s(std::uint16_t key, std::span<const double> v);
    void i64s(std::uint16_t key, std::span<const std::int64_t> v);
    void record(std::uint16_t key, const Writer& nested);

    [[nodiscard]] const std::vector<std::byte>& bytes() const noexcept { return bytes_; }

private:
    void field(std::uint16_t key, std::span<const std::byte> payload);
    std::vector<std::byte> bytes_;
};

// A parsed record. It refers to the bytes it was parsed from, which must outlive it.
class Reader {
public:
    // invalid_input ("corrupt ...") when a field runs past the end of the bytes.
    [[nodiscard]] static std::expected<Reader, base::Error> parse(std::span<const std::byte> bytes);

    [[nodiscard]] bool has(std::uint16_t key) const noexcept;
    // Each getter returns invalid_input when the field is absent or its length is wrong for its
    // type; the *_or forms return the default when it is absent.
    [[nodiscard]] std::expected<std::uint64_t, base::Error> u64(std::uint16_t key) const;
    [[nodiscard]] std::expected<std::int64_t, base::Error> i64(std::uint16_t key) const;
    [[nodiscard]] std::expected<double, base::Error> f64(std::uint16_t key) const;
    [[nodiscard]] std::expected<std::string, base::Error> text(std::uint16_t key) const;
    [[nodiscard]] std::expected<std::vector<double>, base::Error> f64s(std::uint16_t key) const;
    [[nodiscard]] std::expected<std::vector<std::int64_t>, base::Error> i64s(
        std::uint16_t key) const;
    [[nodiscard]] std::expected<Reader, base::Error> record(std::uint16_t key) const;
    // Every field with this key, in order (a list).
    [[nodiscard]] std::expected<std::vector<Reader>, base::Error> records(std::uint16_t key) const;

    [[nodiscard]] std::expected<std::uint64_t, base::Error> u64_or(std::uint16_t key,
                                                                   std::uint64_t fallback) const;
    [[nodiscard]] std::expected<std::int64_t, base::Error> i64_or(std::uint16_t key,
                                                                  std::int64_t fallback) const;
    [[nodiscard]] std::expected<std::vector<double>, base::Error> f64s_or_empty(
        std::uint16_t key) const;

private:
    struct Field {
        std::uint16_t key;
        std::span<const std::byte> bytes;
    };
    [[nodiscard]] std::optional<std::span<const std::byte>> find(std::uint16_t key) const noexcept;
    [[nodiscard]] std::expected<std::span<const std::byte>, base::Error> need(
        std::uint16_t key) const;
    std::vector<Field> fields_;
};

base::Error corrupt(std::string what);

}  // namespace NiTCAD::io::detail
